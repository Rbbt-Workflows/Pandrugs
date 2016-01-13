require 'rbbt-util'
require 'rbbt/resource'

module Pandrugs
  extend Resource
  self.subdir = 'share/databases/Pandrugs'

  def self._pandrug_file
    Rbbt.data["Pandrugs.tsv"].find(:lib)
  end

  Pandrugs.claim Pandrugs.gene_drugs, :proc do 
    fields = %w(source standard_drug_name family status pathology cancer extra extra2 target_marker resistance alteration score)
    tsv = Pandrugs._pandrug_file.tsv :header_hash => "", :merge => true, :zipped => true, :fix => Proc.new{|l| l.gsub(/\s+\|\s+/,',')}, :fields => fields
    tsv.key_field = "Associated Gene Name"
    tsv.to_s
  end

  Pandrugs.claim Pandrugs.drug_ids, :proc do 
    fields = %w(source_drug_name show_drug_name)
    tsv = Pandrugs._pandrug_file.tsv :header_hash => "", :key_field => "standard_drug_name", :fix => Proc.new{|l| l.gsub(/\s+\|\s+/,',')}, :fields => fields
    tsv.to_s
  end
end
