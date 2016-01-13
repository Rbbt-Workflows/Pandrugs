require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

require 'rbbt/sources/organism'
require 'rbbt/sources/Pandrugs'
require 'rbbt/knowledge_base/Pandrugs'

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

end

#require 'Pandrugs/tasks/basic.rb'

#require 'rbbt/knowledge_base/Pandrugs'
#require 'rbbt/entity/Pandrugs'

