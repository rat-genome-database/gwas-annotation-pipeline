package edu.mcw.rgd.gwasAnnot;

public class dbSnp {
    private String snpName;
    private String source;
    private Integer mapKey;
    private String chromosome;
    private Integer position;
    private String snpClass;
    private String molType;
    private String genotype;
    private String hetType;
    private Integer avgHetroScore;
    private Double stdErr;
    private String allele;
    private Double mafFrequency;
    private Integer mafSampleSized;
    private Integer mapLocCount;
    private String mafAllele;
    private Integer dbSnpId;
    private String clinicalSignificance;
    private String refAllele;

    public String getSnpName() {
        return snpName;
    }

    public void setSnpName(String snpName) {
        this.snpName = snpName;
    }

    public String getSource() {
        return source;
    }

    public void setSource(String source) {
        this.source = source;
    }

    public Integer getMapKey() {
        return mapKey;
    }

    public void setMapKey(Integer mapKey) {
        this.mapKey = mapKey;
    }

    public String getChromosome() {
        return chromosome;
    }

    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }

    public Integer getPosition() {
        return position;
    }

    public void setPosition(Integer position) {
        this.position = position;
    }

    public String getSnpClass() {
        return snpClass;
    }

    public void setSnpClass(String snpClass) {
        this.snpClass = snpClass;
    }

    public String getMolType() {
        return molType;
    }

    public void setMolType(String molType) {
        this.molType = molType;
    }

    public String getGenotype() {
        return genotype;
    }

    public void setGenotype(String genotype) {
        this.genotype = genotype;
    }

    public String getHetType() {
        return hetType;
    }

    public void setHetType(String hetType) {
        this.hetType = hetType;
    }

    public Integer getAvgHetroScore() {
        return avgHetroScore;
    }

    public void setAvgHetroScore(Integer avgHetroScore) {
        this.avgHetroScore = avgHetroScore;
    }

    public Double getStdErr() {
        return stdErr;
    }

    public void setStdErr(Double stdErr) {
        this.stdErr = stdErr;
    }

    public String getAllele() {
        return allele;
    }

    public void setAllele(String allele) {
        this.allele = allele;
    }

    public Double getMafFrequency() {
        return mafFrequency;
    }

    public void setMafFrequency(Double mafFrequency) {
        this.mafFrequency = mafFrequency;
    }

    public Integer getMafSampleSized() {
        return mafSampleSized;
    }

    public void setMafSampleSized(Integer mafSampleSized) {
        this.mafSampleSized = mafSampleSized;
    }

    public Integer getMapLocCount() {
        return mapLocCount;
    }

    public void setMapLocCount(Integer mapLocCount) {
        this.mapLocCount = mapLocCount;
    }

    public String getMafAllele() {
        return mafAllele;
    }

    public void setMafAllele(String mafAllele) {
        this.mafAllele = mafAllele;
    }

    public Integer getDbSnpId() {
        return dbSnpId;
    }

    public void setDbSnpId(Integer dbSnpId) {
        this.dbSnpId = dbSnpId;
    }

    public String getClinicalSignificance() {
        return clinicalSignificance;
    }

    public void setClinicalSignificance(String clinicalSignificance) {
        this.clinicalSignificance = clinicalSignificance;
    }

    public String getRefAllele() {
        return refAllele;
    }

    public void setRefAllele(String refAllele) {
        this.refAllele = refAllele;
    }
}
