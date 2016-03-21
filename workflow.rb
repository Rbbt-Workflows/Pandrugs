require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

require 'rbbt/sources/organism'
require 'rbbt/sources/Pandrugs'
require 'rbbt/sources/clinvar'
require 'rbbt/knowledge_base/Pandrugs'

Workflow.require_workflow "Sample"
Workflow.require_workflow "Study"
Workflow.require_workflow "InterPro"
module Pandrugs
  extend Workflow

  input :genes, :array, "List of genes"
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  task :annotate_genes => :tsv do |genes, organism|
    gname = Organism.identifiers(organism).index(:target => "Associated Gene Name", :persist => true).chunked_values_at(genes).compact

    tsv = Pandrugs.gene_drugs.tsv(:persist => true, :persist_update => true, :merge => true,).
            select(:key => gname).
            reorder("standard_drug_name", nil, :zipped => true)

    tsv.namespace = organism
    tsv
  end

  task :gene_cancer_role => :tsv do
    types = TSV.setup({}, :key_field => "Associated Gene Name", :fields => ["Oncogene/Tumor Suppressor"], :type => :single)

    TSV.traverse Rbbt.data["CGC.tsv"].find(:lib), :fields => ["Molecular Genetics"], :type => :single, :into => types, :header_hash => "" do |gene,type|
      n_type = case type 
               when "Dom" 
                 "Oncogene" 
               when "Rec" 
                 "Tumor suppressor"
               else
                 next
               end
      [gene, n_type]
    end

    TSV.traverse Rbbt.data["HCD.oncodriveROLE.0.3-0.7.txt"].find(:lib), :fields => ["oncodriveROLE"], :type => :single, :into => types, :header_hash => "" do |gene,type|
      n_type = case type 
               when "Activating" 
                 "Oncogene" 
               when "Loss of function" 
                 "Tumor suppressor"
               else
                 next
               end
      [gene, n_type]
    end

    types.to_s
  end

  task :COSMIC_incidence => :tsv do
    dumper = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => ["Sample"], :type => :flat
    dumper.init
    TSV.traverse COSMIC.mutations, :key_field => "Genomic Mutation", :fields => ["Sample name"], :bar => true, :into => dumper do |mut,values|
      [mut, values]
    end
    TSV.collapse_stream dumper.stream
  end

  dep :COSMIC_incidence
  task :COSMIC_gene_incidence => :tsv do
    organism = COSMIC.organism
    Step.wait_for_jobs dependencies

    io = CMD.cmd("cut -f 1 '#{step(:COSMIC_incidence).path}'", :pipe => true)
    job = Sequence.job(:mutated_isoforms_fast, nil, :mutations => io, :principal => true, :non_synonymous => true, :organism => organism)
    job.run(true)

    incidence = step(:COSMIC_incidence).load

    index = Organism.identifiers(organism).index :persist => true

    dumper = TSV::Dumper.new :key_field => "Ensembl Gene ID", :fields => ["Sample"], :type => :flat
    dumper.init
    TSV.traverse job, :into => dumper, :bar => true do |mut, mis|
      mut = mut.first if Array === mut
      samples = incidence[mut]

      genes = mis.collect do |mi|
        index[mi.partition(":").first]
      end.uniq.compact

      res = genes.collect do |gene|
        [gene, samples]
      end

      res.extend MultipleResult
      res
    end
    TSV.collapse_stream dumper.stream
  end

  dep :COSMIC_gene_incidence
  dep :COSMIC_incidence
  task :COSMIC_counts => :tsv do
    Misc.open_pipe do |sin|
      dumper_mut = TSV::Dumper.new :key_field => "Thing", :fields => ["Sample count"], :type => :single
      dumper_mut.init
      max = 0
      TSV.traverse step(:COSMIC_incidence), :into => dumper_mut, :bar => true do |mut,samples|
        mut = mut.first if Array === mut
        count = samples.length
        max = count if count > max
        [mut, count]
      end
      Misc.consume_stream(dumper_mut.stream,false, sin)
      set_info :max_mutations, max

      dumper_gene = TSV::Dumper.new :key_field => "Thing", :fields => ["Sample count"], :type => :single
      max = 0
      TSV.traverse step(:COSMIC_gene_incidence), :into => dumper_gene, :bar => true do |gene,samples|
        gene = gene.first if Array === gene
        count = samples.length
        max = count if count > max
        [gene, count]
      end
      Misc.consume_stream(dumper_gene.stream,false, sin)
      set_info :max_genes, max
    end
  end

  dep Sample, :sequence_ontology, :principal => true, :non_synonymous => true
  dep Sample, :DbNSFP
  dep Sample, :genomic_mutation_consequence, :principal => true, :non_synonymous => true
  task :mi_info4score => :tsv do
    Step.wait_for_jobs dependencies
    tsv = step(:sequence_ontology).load.to_double
    tsv = tsv.attach(step(:genomic_mutation_consequence).load.to_flat)
    tsv = tsv.attach(step(:DbNSFP))
    organism = step(:organism).load
    index = Organism.identifiers(organism).index :persist => true
    tsv.add_field "Ensembl Gene ID" do |mi, values|
      next unless mi =~ /ENSP/
      index[mi.partition(":").first]
    end
  end

  dep :mi_info4score
  dep Sample, :homozygous
  dep Sample, :mi
  dep Sample, :organism
  dep InterPro, :domains, :mutated_isoforms => :mi, :organism => :organism
  dep Sample, :genomic_mutation_annotations
  dep Sample, :mutation_info
  dep :gene_cancer_role do |jobname,options|
    Pandrugs.job(:gene_cancer_role, nil)
  end
  task :gene_score => :tsv do
    
    organism = step(:organism).load

    clinvar = ClinVar.mi_summary.tsv :fields => ["ClinicalSignificance"], :type => :single, :persist => true
    essenciality = Rbbt.data["gene_essentiality_score.tsv"].tsv :header_hash => '', :fields => ["max_min_score"], :type => :single, :cast => :to_f
    name_index = Organism.identifiers(organism).index :persist => true, :target => "Associated Gene Name"
    ensp2ensg = Organism.identifiers(organism).index :persist => true, :target => "Ensembl Gene ID", :fields => ["Ensembl Protein ID"]

    gene_cancer_role = step(:gene_cancer_role).load

    domains = step(:domains).join.path.tsv :fields => ["InterPro Domain"], :type => :flat
    intrp2pfam = InterPro.pfam_equivalences.index :target => "Pfam Domain", :persist => true

    good_pfam_domains = Rbbt.data["domains.tsv"].tsv(:key_field => ["Domain code"], :fields => [], :type => :list).keys

    homozygous = step(:homozygous).load

    count_job = Pandrugs.job(:COSMIC_counts)
    count_job.produce
    cosmic_counts = count_job.path.tsv :persist => true, :persist_file => count_job.file(:persistence)
    max_muts = count_job.info[:max_mutations]
    max_genes = count_job.info[:max_genes]

    organism = step(:organism).load
    dumper = TSV::Dumper.new :key_field => "Ensembl Gene ID", :fields => ["Gene Score"], :type => :single, :namespace => organism
    dumper.init
    pasted = TSV.paste_streams([step(:mi_info4score), step(:genomic_mutation_annotations), step(:mutation_info)])
    parser = TSV::Parser.new pasted, :type => :double
    fields = parser.fields
    TSV.traverse parser, :into => dumper, :bar => true do |mutation, values|
      values = values.collect{|v| v.compact.reject{|s| s.empty?}.first }
      mutation = mutation.first if Array === mutation

      info = Hash[*fields.zip(values).flatten]

      mi = info["Mutated Isoform"]
      if mi
        protein = mi.partition(":").first
        gene = ensp2ensg[protein]
      else
        gene = info["Ensembl Gene ID"]
      end
      gene_name = name_index[gene]

      role = gene_cancer_role[gene_name] || "Oncogene"
      oncogene = role == 'Oncogene'

      score = 0

      # Damage prediction scores
      score += 0.125/3 if info["Polyphen2_HDIV_score"].to_f > 0.435

      score += 0.125/3 if info["SIFT_score"].to_f > 0.435

      score += 0.125/3 if info["MetaSVM_score"].to_f > 0  #instead of condel

      score += (oncogene ? 0.125/3 : 0.03125) if info["FATHMM_score"].to_f < -1.5 and info["FATHMM_score"] != -999


      # COSMIC Freq
      mut_samples = cosmic_counts[mutation]
      gene_samples = cosmic_counts[gene]

      if mut_samples and oncogene
        mut_samples = mut_samples.to_i
        if mut_samples > 100
          score += 0.125/3
        else
          score += 0.125/3 * Math.log(mut_samples)/Math.log(max_muts)
        end
      end

      if gene_samples
        gene_samples = gene_samples.to_i
        if gene_samples > 100
          score += oncogene ? 0.125/3 : 0.03125
        else
          score += (oncogene ? 0.125/3 : 0.03125) * Math.log(gene_samples)/Math.log(max_genes)
        end
      end

      # MI type scores
      score += 0.125 if %w(stop_gain frameshift_variant missense_variant inframe_insertion inframe_deletion).include?(info["SO Term"])

      #ExAC
      score += 0.125/2 if info["ExAC_Adj_AF"].to_f < 0.01 

      #GMAF
      score += 0.125/2 if info["1000Gp3_AF"].to_f < 0.01 

      # Clinvar
      score += 0.125 if clinvar[mi] == "Pathogenic"

      # Zygosity
      score += (oncogene ? 0.125 : 0.1875) if homozygous.include? mutation

      # Esenciality
      score += 0.125 * (essenciality[gene_name] || 0)

      # Domains
      ip_domains = domains[mi] || []
      pfam_domains = intrp2pfam.values_at(*ip_domains).compact
      
      if pfam_domains.any?
        if (good_pfam_domains & pfam_domains).any?
          score += 0.125 
        else
          score += 0.125/2
        end
      end

      [gene, score]
    end

    TSV.open(dumper.stream, :type => :double, :merge => true).to_list do |values|
      values.max
    end
  end

end

module Sample
  dep :gene_mutation_status
  dep Pandrugs, :gene_score
  task :pandrugs => :tsv do

    gene_status = step(:gene_mutation_status).load
    affected_genes = gene_status.select("affected" => "true").keys
    broken_genes = gene_status.select("broken" => "true").keys

    broken_genes = Translation.translate(organism, "Associated Gene Name", broken_genes)

    tsv = Pandrugs.job(:annotate_genes, sample, :genes => affected_genes, :organism => organism).run

    tsv.add_field "Broken" do |drug,values|
      broken_genes.include? values.first
    end
    scores = step(:gene_score).load

    scores.identifiers = Organism.identifiers(organism)
    scores.change_key("Associated Gene Name")
    tsv = tsv.attach scores, :identifiers => Organism.identifiers(organism)

    tsv
  end
end


module Study

  extend Workflow

  dep Sample, :pandrugs do |jobname,options|
    study = Study.setup(jobname.dup)
    jobs = study.genotyped_samples.collect{|sample| Sample.setup(sample, :cohort => study); sample.pandrugs(:job, options) }.flatten
    Misc.bootstrap(jobs, 3, :bar => "Processing gene_sample_mutation_status", :respawn => :always) do |job|
      job.produce
      nil
    end
    jobs
  end

  task :pandrugs => :tsv do
    Step.wait_for_jobs dependencies
    parser = TSV::Parser.new dependencies.first
    fields = parser.fields
    fields.unshift "Sample"
    header = TSV.header_lines(parser.key_field, parser.fields, parser.options.merge(:type => :double))

    io = Misc.open_pipe do |sin|
      sin.puts header

      TSV.traverse dependencies, :type => :array do |job|
        sample = job.clean_name.split(":").last
        TSV.traverse job, :type => :array do |line|
          next if line =~ /^#/
            gene,*rest = line.split("\t")
          parts = [gene, sample]
          parts.concat rest
          sin.puts parts * "\t"
        end
      end
    end

    TSV.collapse_stream io
  end
end

Study.update_task_properties
Sample.update_tasks_property_bindings

#
##require 'Pandrugs/tasks/basic.rb'
#
##require 'rbbt/knowledge_base/Pandrugs'
##require 'rbbt/entity/Pandrugs'
#

