package edu.mcw.rgd;

import edu.mcw.rgd.dao.impl.*;
import edu.mcw.rgd.dao.impl.variants.VariantDAO;
import edu.mcw.rgd.dao.spring.StringMapQuery;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.datamodel.ontology.Annotation;
import edu.mcw.rgd.datamodel.ontologyx.Term;
import edu.mcw.rgd.datamodel.ontologyx.TermSynonym;
import edu.mcw.rgd.datamodel.variants.VariantMapData;
import edu.mcw.rgd.process.Utils;

import java.util.*;

public class DAO {

    AnnotationDAO adao = new AnnotationDAO();
    GWASCatalogDAO gdao = new GWASCatalogDAO();
    OntologyXDAO odao = new OntologyXDAO();
    VariantDAO vdao = new VariantDAO();
    QTLDAO qdao = new QTLDAO();
    NotesDAO noteDAO = new NotesDAO();
    RGDManagementDAO managementDAO = new RGDManagementDAO();

    public String getConnectionInfo() {
        return adao.getConnectionInfo();
    }

    public List<Annotation> getAnnotations(int annotatedObjectRGDId, String termAcc, int createdBy) throws Exception {
        return adao.getAnnotations(annotatedObjectRGDId, termAcc, createdBy);
    }

    public List<Annotation> getAnnotationsModifiedBeforeTimestampForEFO(Date dt, int createdBy) throws Exception{
        return adao.getAnnotationsModifiedBeforeTimestamp(createdBy, dt, "T");
    }

    public List<Annotation> getAnnotationsModifiedBeforeTimestampForCMO(Date dt, int createdBy) throws Exception{
        return adao.getAnnotationsModifiedBeforeTimestamp(createdBy, dt, "L");
    }

    public List<Annotation> getAnnotationsModifiedBeforeTimestampForVT(Date dt, int createdBy) throws Exception{
        return adao.getAnnotationsModifiedBeforeTimestamp(createdBy, dt, "V");
    }

    public List<Annotation> getAnnotationsModifiedBeforeTimestampForRDO(Date dt, int createdBy) throws Exception{
        return adao.getAnnotationsModifiedBeforeTimestamp(createdBy,dt,"D");
    }
    public List<Annotation> getAnnotationsModifiedBeforeTimestamp(Date dt, int createdBy) throws Exception {
        List<Annotation> annotations = new ArrayList<>();
        List<Annotation> efo = getAnnotationsModifiedBeforeTimestampForEFO(dt, createdBy);
        List<Annotation> cmo = getAnnotationsModifiedBeforeTimestampForCMO(dt,createdBy);
        List<Annotation> vt = getAnnotationsModifiedBeforeTimestampForVT(dt, createdBy);
        if (efo!=null && !efo.isEmpty())
            annotations.addAll(efo);
        if (cmo!=null && !cmo.isEmpty())
            annotations.addAll(cmo);
        if (vt!=null && !vt.isEmpty())
            annotations.addAll(vt);
        return annotations;
    }

    public List<GWASCatalog> getGWASCatalog() throws Exception{
        return gdao.getFullCatalog();
    }

    public Term getTermByName(String term, String ontId) throws Exception{
        return odao.getTermByTermName(term,ontId);
    }

    public Term getTermByAccId(String accId) throws Exception{
        return odao.getTermByAccId(accId);
    }

    public int GenerateNextQTLSeqForGwas()throws Exception{
        return adao.getNextKeyFromSequence("GWAS_QTL_SEQ");
    }

    public Integer getVariantRgdId(GWASCatalog gc) throws Exception{
        List<VariantMapData> vmds = vdao.getAllVariantByRsIdAndMapKey(gc.getSnps(),38);
        for (VariantMapData v : vmds){
            if (v.getVariantNucleotide().equals(gc.getStrongSnpRiskallele()))
                return (int)v.getId();
        }
        return null;
    }

    public List<Integer> getVariantRgdIds(GWASCatalog gc) throws Exception{
        List<VariantMapData> vmds = vdao.getAllVariantByRsIdAndMapKey(gc.getSnps(),38);
        List<Integer> rgdIds = new ArrayList<>();
        for (VariantMapData v : vmds){
            if (!rgdIds.contains((int) v.getId()))
                rgdIds.add((int)v.getId());
        }
        return rgdIds;
    }

    public List<TermSynonym> getTermSynonymsBySynonymName(String synonym) throws Exception {
        return odao.getTermSynonymBySynonymName(synonym);
    }

    public QTL getQtlByRgdId(int rgdId) throws Exception{
        return qdao.getQTL(rgdId);
    }

    public void insertQTL(QTL q) throws Exception {
        qdao.insertQTL(q,"ACTIVE",1);
    }

    public int insertAnnotation(Annotation a) throws Exception {
        return adao.insertAnnotation(a);
    }

    public List<Term> getQTLTraitTerms(int rgdId) throws Exception{
        List<Note> notes = new ArrayList<>();
        List<Term> traitTerms = getTermsForObject(rgdId, "V");
        if( traitTerms.isEmpty() ) {
            notes = noteDAO.getNotes(rgdId, "qtl_trait");
            if( !notes.isEmpty() ) {
                Term traitTerm = new Term();
                traitTerm.setTerm(notes.get(0).getNotes());
                traitTerms.add(traitTerm);
            }
        }
        return traitTerms;
    }

    public List<Term> getQTLMeasurementTerms(QTL q) throws Exception{
        List<Note> notes = new ArrayList<>();
        List<Term> measurementTerms = getTermsForObject(q.getRgdId(),"L");
        if( measurementTerms.isEmpty() ) {
            notes = noteDAO.getNotes(q.getRgdId(), "qtl_subtrait");
            if( !notes.isEmpty() ) {
                Term measurementTerm = new Term();
                measurementTerm.setTerm(notes.get(0).getNotes());
                measurementTerms.add(measurementTerm);
            }
        } else if( measurementTerms.size()>1 ) {
            // multiple CMO terms -- pick the most significant CMO term if available
            if( !Utils.isStringEmpty(q.getMostSignificantCmoTerm()) ) {
                Term term = odao.getTermWithStatsCached(q.getMostSignificantCmoTerm());
                if( term!=null ) {
                    measurementTerms.clear();
                    measurementTerms.add(term);
                }
            }
        }
        return measurementTerms;
    }

    private List<Term> getTermsForObject(int rgdId, String aspect) throws Exception {
        List<Term> terms = new ArrayList<>();
        for( StringMapQuery.MapPair pair: adao.getAnnotationTermAccIds(rgdId, aspect) ) {
            Term term = new Term();
            term.setAccId(pair.keyValue);
            term.setTerm(pair.stringValue);
            terms.add(term);
        }
        if( terms.size()>1 ) {
            Collections.sort(terms, new Comparator<Term>() {
                @Override
                public int compare(Term o1, Term o2) {
                    return Utils.stringsCompareToIgnoreCase(o1.getTerm(), o2.getTerm());
                }
            });
        }
        return terms;
    }

    public RgdId createRgdId(int objectKey, String objectStatus, String notes, int mapKey) throws Exception{
        int speciesKey = SpeciesType.getSpeciesTypeKeyForMap(mapKey);
        return managementDAO.createRgdId(objectKey, objectStatus, notes, speciesKey);
    }

    public int updateGwasQtlRgdIdBatch(Collection<GWASCatalog> update) throws Exception{
        return gdao.updateGwasQtlRgdIdBatch(update);
    }
}