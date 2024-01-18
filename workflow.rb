require 'rbbt-util'
require 'rbbt/workflow'


Misc.add_libdir if __FILE__ == $0

require 'rbbt/sources/organism'
require 'rbbt/sources/Pandrugs'
require 'rbbt/sources/clinvar'
require 'rbbt/sources/pfam'
require 'rbbt/knowledge_base/Pandrugs'

Workflow.require_workflow "Sample"
Workflow.require_workflow "Genomics"
Workflow.require_workflow "Study"
Workflow.require_workflow "InterPro"

module Pandrugs
  extend Workflow

  input :genes, :array, "List of genes"
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  task :annotate_genes => :tsv do |genes, organism|
    gname = Organism.identifiers(organism).index(:target => "Associated Gene Name", :persist => true).chunked_values_at(genes).compact

    matches = Pandrugs.gene_drugs.tsv(:persist => true, :merge => true, :type => :double).select(:key => gname)

    tsv = matches.reorder("standard_drug_name", nil, :zipped => true, :merge => true, :persist => false)

    tsv.fields = tsv.fields.collect{|f| f == 'score' ? 'Drug Score' : f }

    tsv.namespace = organism
    tsv
  end

  task :gene_cancer_role => :tsv do
    types = TSV.setup({}, :key_field => "Associated Gene Name", :fields => ["Oncogene/Tumor Suppressor"], :type => :single)

    TSV.traverse Rbbt.root.data["CGC.tsv"].find(:lib), :fields => ["Molecular Genetics"], :type => :single, :into => types, :header_hash => "" do |gene,type|
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

    TSV.traverse Rbbt.root.data["HCD.oncodriveROLE.0.3-0.7.txt"].find(:lib), :fields => ["oncodriveROLE"], :type => :single, :into => types, :header_hash => "" do |gene,type|
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

  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  task :COSMIC_incidence => :tsv do |organism|
    dumper = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => ["Sample"], :type => :flat
    dumper.init
    build = Organism.GRC_build(organism)
    TSV.traverse COSMIC[build].mutations, :key_field => "Genomic Mutation", :fields => ["Sample name"], :bar => true, :into => dumper do |mut,values|
      [mut, values]
    end
    TSV.collapse_stream dumper.stream
  end

  dep :COSMIC_incidence, :compute => :produce
  task :COSMIC_gene_incidence => :tsv do
    incidence = step(:COSMIC_incidence).load
    organism = step(:COSMIC_incidence).inputs[:organism]

    io = CMD.cmd("cut -f 1 '#{step(:COSMIC_incidence).path}'", :pipe => true)
    job = Sequence.job(:mutated_isoforms_fast, nil, :mutations => io, :principal => true, :non_synonymous => true, :organism => organism)
    job.produce


    index = Organism.identifiers(organism).index :persist => true

    dumper = TSV::Dumper.new :key_field => "Ensembl Gene ID", :fields => ["Sample"], :type => :double
    dumper.init
    TSV.traverse job, :into => dumper, :bar => true do |mut, mis|
      mut = mut.first if Array === mut
      samples = incidence[mut]

      genes = mis.collect do |mi|
        index[mi.partition(":").first]
      end.uniq.compact

      samples = samples.collect{|s| s.split("|")}
      res = genes.collect do |gene|
        [gene, samples]
      end

      res.extend MultipleResult
      res
    end
    TSV.collapse_stream dumper.stream
  end

  dep :COSMIC_gene_incidence, :compute => :produce
  task :COSMIC_counts => :tsv do
    Misc.open_pipe do |sin|
      dumper_mut = TSV::Dumper.new :key_field => "Thing", :fields => ["Sample count"], :type => :single
      dumper_mut.init
      max = 0
      cosmic_incidence = step(:COSMIC_gene_incidence).step(:COSMIC_incidence)
      TSV.traverse cosmic_incidence, :into => dumper_mut, :bar => true do |mut,values|
        samples = values.first
        mut = mut.first if Array === mut
        count = samples.length
        max = count if count > max
        [mut, count]
      end
      Misc.consume_stream(dumper_mut.stream,false, sin, false)
      set_info :max_mutations, max

      dumper_gene = TSV::Dumper.new :key_field => "Thing", :fields => ["Sample count"], :type => :single
      max = 0
      TSV.traverse step(:COSMIC_gene_incidence), :into => dumper_gene, :bar => true do |gene,values|
        samples = values.first
        gene = gene.first if Array === gene
        count = samples.length
        max = count if count > max
        [gene, count]
      end
      Misc.consume_stream(dumper_gene.stream,false, sin, false)
      set_info :max_genes, max
    end
  end

  dep Sample, :sequence_ontology, :principal => true, :non_synonymous => true, :compute => :produce
  dep Sample, :DbNSFP, :compute => :produce
  dep Sample, :genomic_mutation_splicing_consequence, :principal => true, :compute => :produce
  task :mi_info4score => :tsv do
    tsv = step(:sequence_ontology).load.to_double

    tsv = tsv.attach(step(:DbNSFP).path)

    spl = step(:genomic_mutation_splicing_consequence)
    tsv = tsv.attach(spl)

    organism = step(:organism).load
    ensp2ensg = Organism.identifiers(organism).index :persist => true
    enst2ensg = Organism.transcripts(organism).index :persist => true, :target => "Ensembl Gene ID", :fields => ["Ensembl Transcript ID"]
    tsv.add_field "Ensembl Gene ID" do |mut, values|
      mis = values["Mutated Isoform"]

      mi_genes = values["Mutated Isoform"].collect{|mi|
        next unless mi =~ /ENSP/
          ensp2ensg[mi.partition(":").first]
      }.compact

      splicing_genes = values["Affected Transcripts"].collect{|t|
        enst2ensg[t]
      }.compact

      (mi_genes + splicing_genes).uniq * "|"
    end
  end

  dep :mi_info4score, :compute => :produce
  dep Sample, :homozygous, :compute => :produce
  dep Sample, :mi, :compute => :produce
  dep Sample, :organism, :compute => :produce
  dep InterPro, :domains, :mutated_isoforms => :mi, :organism => :organism, :compute => :produce
  dep Sample, :annotate_Genomes1000, :compute => :produce
  dep Sample, :annotate_DbSNP, :compute => :produce
  dep Sample, :mutation_info, :compute => :produce
  dep :gene_cancer_role, :compute => :produce do |jobname,options|
    Pandrugs.job(:gene_cancer_role, nil)
  end
  input :use_COSMIC, :boolean, "Use COSMIC to score variants", false
  task :gene_score => :tsv do |use_COSMIC|

    organism = step(:organism).load

    build = Organism.hg_build(organism)

    clinvar = ClinVar[build].mi_summary.tsv :fields => ["ClinicalSignificance"], :type => :single, :persist => true
    essenciality = Rbbt.root.data["gene_essentiality_score.tsv"].tsv :header_hash => '', :fields => ["max_min_score"], :type => :single, :cast => :to_f
    name_index = Organism.identifiers(organism).index :persist => true, :target => "Associated Gene Name"
    ensp2ensg = Organism.identifiers(organism).index :persist => true, :target => "Ensembl Gene ID", :fields => ["Ensembl Protein ID"]

    gene_cancer_role = step(:gene_cancer_role).load

    domains = step(:domains).join.path.tsv :fields => ["InterPro Domain"], :type => :flat
    intrp2pfam = InterPro.pfam_equivalences.index :target => "Pfam Domain", :persist => true

    good_pfam_domains = Rbbt.root.data["domains.tsv"].tsv(:key_field => ["Domain code"], :fields => [], :type => :list).keys

    homozygous = step(:homozygous).load

    if use_COSMIC
      count_job = Pandrugs.job(:COSMIC_counts)
      count_job.produce
      cosmic_counts = count_job.path.tsv :persist => true, :persist_file => count_job.file(:persistence)
      max_muts = count_job.info[:max_mutations]
      max_genes = count_job.info[:max_genes]
    end

    organism = step(:organism).load
    dumper = TSV::Dumper.new :key_field => "Ensembl Gene ID", :fields => ["Gene Score"], :type => :single, :namespace => organism
    dumper.init
    pasted = TSV.paste_streams([step(:mi_info4score), step(:annotate_Genomes1000), step(:annotate_DbSNP), step(:mutation_info)])
    parser = TSV::Parser.new pasted, :type => :double
    fields = parser.fields
    gene_pos = fields.index "Ensembl Gene ID"
    TSV.traverse parser, :into => dumper, :bar => true do |mutation, lvalues|
      mutation = mutation.first if Array === mutation
      lvalues = [lvalues[gene_pos]] + lvalues

      res = Misc.zip_fields(lvalues).collect do |values|
        #values = values.collect{|v| v.compact.reject{|s| s.empty?}.first }

        info = {}
        fields.zip(values[1..-1]).each do |k,v|
          info[k] = v if info[k].nil? or info[k].empty?
        end

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
        if use_COSMIC
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
      res.extend MultipleResult
      res
    end

    TSV.open(dumper.stream, :type => :double, :merge => true).to_list do |values|
      values.max
    end
  end

end

module Sample

  dep :organism
  dep :affected_genes, :principal => true, :compute => :produce
  dep Pandrugs, :annotate_genes, :genes => :affected_genes, :organism => :organism, :compute => :produce
  task :pandrugs_annotated_genes => :tsv do
    TSV.get_stream step(:annotate_genes)
  end

  dep :pandrugs_annotated_genes
  dep :gene_mutation_status, :principal => true, :compute => :produce
  dep Pandrugs, :gene_score
  task :pandrugs => :tsv do

    gene_status = step(:gene_mutation_status).load
    affected_genes = gene_status.select("affected" => "true").keys
    tsv = step(:pandrugs_annotated_genes).load
    tsv.swap_id("Associated Gene Name", "Ensembl Gene ID", :persist => true)

    scores = step(:gene_score).load

    scores.identifiers = Organism.identifiers(organism)
    tsv = tsv.attach scores, :identifiers => Organism.identifiers(organism), :persist_input => true, :one2one => true

    tsv
  end
end


module Study

  dep Sample, :pandrugs, :compute => :bootstrap do |jobname,options|
    study = Study.setup(jobname.dup)
    study.genotyped_samples.collect{|sample| Sample.setup(sample, :cohort => study); sample.pandrugs(:job, options) }.flatten
  end
  task :pandrugs => :tsv do
    parser = TSV::Parser.new dependencies.first
    fields = parser.fields
    fields.unshift "Sample"
    header = TSV.header_lines(parser.key_field, parser.fields, parser.options.merge(:type => :double))

    io = Misc.open_pipe do |sin|
      sin.puts header

      TSV.traverse dependencies, :type => :array, :bar => "Joining sample information" do |job|
        sample = job.clean_name.split(":").last
        TSV.traverse job, :type => :array do |line|
          next if line =~ /^#/

          gene,*rest = line.split("\t")
          rest_split = rest.collect{|r| r.empty? ? [" "] : r.split("|",-1).collect{|v| v.empty? ? " " : v}}

          if rest_split.collect{|p| p.length}.uniq.length > 1
            raise "Number of fields does not match in:\nline"
          end

          num = rest_split.first.length
          sample_str = ([sample] * num) * "|"

          parts = [gene, sample_str]
          parts.concat rest_split.collect{|v| v*"|"}

          sin.puts parts * "\t"
        end
      end
    end

    TSV.collapse_stream io
  end
end

Study.update_task_properties
Sample.update_task_properties

#
##require 'Pandrugs/tasks/basic.rb'
#
##require 'rbbt/knowledge_base/Pandrugs'
##require 'rbbt/entity/Pandrugs'
#

