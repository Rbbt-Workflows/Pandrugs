require 'rbbt/knowledge_base'
require 'rbbt/sources/organism'

$LOAD_PATH.unshift(Rbbt.lib.find(:lib))
require 'rbbt/sources/Pandrugs'

module Pandrugs

  class << self 
    attr_accessor :knowledge_base_dir
  end
  self.knowledge_base_dir = Rbbt.var.knowledge_base.Pandrugs

  def self.organism
    Organism.default_code("Hsa")
  end

  def self.knowledge_base
    @knowledge_base ||= begin
                          kb = KnowledgeBase.new self.knowledge_base_dir, self.organism
                          kb.register "gene_drugs", Pandrugs.gene_drugs, :source => "Associated Gene Name", :target => "standard_drug_name", :merge => true
                          kb
                        end
  end
end
