package edu.mcw.rgd.gwasAnnot;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.datamodel.ontology.Annotation;
import edu.mcw.rgd.datamodel.ontologyx.Term;
import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;

public class RatAnnot {

    private String version;
    private int createdBy;
    private int refRgdId;
    private int refKey;
    private String deleteThresholdForStaleAnnotations;

    Logger status = LogManager.getLogger("ratStatus");
    Logger logStatus = LogManager.getLogger("deleteAnnots");
    SimpleDateFormat sdt = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
    DAO dao = new DAO();
    public void run() throws Exception{
        long time0 = System.currentTimeMillis();
        Date date0 = new Date();

        status.info("   "+dao.getConnectionInfo());

        status.info("GWAS Annotation Pipeline started at "+sdt.format(date0));
        List<GWASCatalog> gwas = dao.getGWASByMapKey(372);
        try {
            runAnnot(gwas);
        }
        catch (Exception e){
            Utils.printStackTrace(e,status);
        }
        status.info("\nTotal pipeline runtime -- elapsed time: "+
                Utils.formatElapsedTime(time0,System.currentTimeMillis()));
    }

    void runAnnot(List<GWASCatalog> gwas) throws Exception {
        HashMap<String, QTL> qtlHashMap = new HashMap<>(); // rsId + PVal is key to help prevent creating duplicates
        List<QTL> existingQtl = new ArrayList<>();
        List<QTL> newQtls = new ArrayList<>();
        List<Annotation> allAnnots = new ArrayList<>();
        List<XdbId> newXdbs = new ArrayList<>();
        List<Integer> qtlRgdIds = new ArrayList<>();
        List<GWASCatalog> update = new ArrayList<>();
        HashMap<Integer, List<String>> rgdToTerm = new HashMap<>();
        HashMap<String, Term> termMap = new HashMap<>();
        HashMap<String, Integer> duplicateCatcher = new HashMap<>();

        for (GWASCatalog g : gwas){
            QTL gwasQtl = new QTL();
            if (Utils.isStringEmpty(g.getSnps()))
                continue;
            if (Utils.isStringEmpty(g.getEfoId()))
                continue;

            if (qtlHashMap.get(g.getVariantRgdId()+"|"+g.getpVal())!=null){
                gwasQtl = qtlHashMap.get(g.getSnps()+"|"+g.getpVal());
            }
            else if (g.getQtlRgdId() != null  && !g.getQtlRgdId().equals(0)) {
                gwasQtl = dao.getQtlByRgdId(g.getQtlRgdId());
                qtlHashMap.put(g.getSnps()+"|"+g.getpVal(), gwasQtl);
                existingQtl.add(gwasQtl);
            }
            else {
                int qtlNum = dao.GenerateNextQTLSeqForRatGwas();
                gwasQtl.setSymbol("GWAS"+qtlNum+"_R");
                gwasQtl.setName(g.getMapTrait() + " QTL GWAS"+ qtlNum);
                gwasQtl.setChromosome(g.getChr());
                RgdId r = dao.createRgdId(RgdId.OBJECT_KEY_QTLS, "ACTIVE", "created by GWAS annotation Pipeline", g.getMapKey());
                gwasQtl.setRgdId(r.getRgdId());
//                gwasQtl.setRgdId(0);
                gwasQtl.setPeakRsId(g.getSnps());
                Double pVal = g.getpVal().doubleValue();
                gwasQtl.setPValue(pVal);
                g.setQtlRgdId(r.getRgdId());
                update.add(g);
                newQtls.add(gwasQtl);
                qtlRgdIds.add(r.getRgdId());
                qtlHashMap.put(g.getVariantRgdId()+"|"+g.getpVal(), gwasQtl);
            }
            if (g.getQtlRgdId()==null || g.getQtlRgdId()==0){
                g.setQtlRgdId(gwasQtl.getRgdId());
                update.add(g);
            }
//            List<XdbId> xdbs = dao.getGwasXdbs(gwasQtl.getRgdId());
//            if (xdbs.isEmpty()){
//                XdbId x = createXdb(g, gwasQtl);
//                newXdbs.add(x);
//            }
            String[] terms = g.getEfoId().split(",");
            rgdToTerm.computeIfAbsent(gwasQtl.getRgdId(), k -> new ArrayList<>());
            rgdToTerm.computeIfAbsent(g.getVariantRgdId(), k -> new ArrayList<>());
            List<String> termsQtl = rgdToTerm.get(gwasQtl.getRgdId());
            List<String> termsVar = rgdToTerm.get(g.getVariantRgdId());


            for (String term : terms){
                Term t = new Term();
                if ( (t = termMap.get(term)) == null ) {
                    t = dao.getTermByAccId(term);
                    termMap.put(term,t);
                }
                Annotation a = new Annotation();
                if (t != null && !termsQtl.contains(t.getAccId()) && !checkAnnotationExist(gwasQtl.getRgdId(), t)){

                    a.setCreatedBy(getCreatedBy());
                    a.setLastModifiedBy(getCreatedBy());
                    a.setAnnotatedObjectRgdId(gwasQtl.getRgdId());
                    a.setRefRgdId(refRgdId);
                    if (t.getAccId().startsWith("CMO"))
                        a.setAspect("L");
                    else // it is VT
                        a.setAspect("V");
                    a.setCreatedDate(new Date());
                    a.setLastModifiedDate(a.getCreatedDate());
                    a.setTerm(t.getTerm());
                    a.setTermAcc(t.getAccId());
                    a.setWithInfo(gwasQtl.getPeakRsId());
                    a.setObjectSymbol(gwasQtl.getSymbol());
                    a.setObjectName(gwasQtl.getName());
                    a.setSpeciesTypeKey(3);
                    a.setDataSrc("RAT_GWAS");
                    a.setEvidence("IAGP");
                    a.setRgdObjectKey(6); // 6 - qtls
                    a.setXrefSource("P50DA037844");
                    String annot = a.getRefRgdId()+"|"+a.getAnnotatedObjectRgdId()+"|"+a.getTermAcc()+"|"+
                            a.getXrefSource()+"|"+a.getQualifier()+"|"+a.getWithInfo()+"|"+a.getEvidence();
                    // copy annot and create for variant
                    if (duplicateCatcher.get(annot)==null) {
                        allAnnots.add(a);
                        duplicateCatcher.put(annot,6);
                    }
//                    terms.add(t.getAccId());
                }
                if (t != null && !termsVar.contains(t.getAccId()) && !checkAnnotationExist(g.getVariantRgdId(), t)) {// same check but for var
                    Annotation aVar = new Annotation();
                    aVar.setAnnotatedObjectRgdId(g.getVariantRgdId());
                    aVar.setCreatedBy(getCreatedBy());
                    aVar.setLastModifiedBy(getCreatedBy());
                    aVar.setRefRgdId(refRgdId);
                    if (t.getAccId().startsWith("CMO"))
                        aVar.setAspect("L");
                    else // it is VT
                        aVar.setAspect("V");
                    aVar.setCreatedDate(new Date());
                    aVar.setTerm(t.getTerm());
                    aVar.setTermAcc(t.getAccId());
                    aVar.setObjectSymbol(g.getSnps());
                    aVar.setSpeciesTypeKey(3);
                    aVar.setDataSrc("RAT_GWAS");
                    aVar.setEvidence("IAGP");
                    aVar.setRgdObjectKey(7);
                    aVar.setXrefSource("P50DA037844");
                    String annot = aVar.getRefRgdId()+"|"+aVar.getAnnotatedObjectRgdId()+"|"+aVar.getTermAcc()+"|"+
                            aVar.getXrefSource()+"|"+aVar.getQualifier()+"|"+aVar.getWithInfo()+"|"+aVar.getEvidence();
                    // copy annot and create for variant
                    if (duplicateCatcher.get(annot)==null) {
                        allAnnots.add(aVar);
                        duplicateCatcher.put(annot,7);
                    }
//                    allAnnots.add(aVar);
                }
            }
            if (!qtlRgdIds.contains(gwasQtl.getRgdId()) && !checkRefAssocExist(gwasQtl.getRgdId())){
                qtlRgdIds.add(gwasQtl.getRgdId());
            }

        } // end gwas for
        if (!existingQtl.isEmpty())
            status.info("\tGWAS QTLs already existing: "+existingQtl.size());
        if (!newQtls.isEmpty()) {
            status.info("\tNew QTLs being made for GWAS: "+newQtls.size());
            dao.insertQTLBatch(newQtls);
        }
        if (!update.isEmpty())
        {
            status.info("\tUpdating GWAS objects: " + update.size());
            dao.updateGwasQtlRgdIdBatch(update);
        }
//        if (!newXdbs.isEmpty())
//        {
//            status.info("\tNew XdbIds being made for QTLs: "+newXdbs.size());
//            dao.insertGwasXdbs(newXdbs);
//        }
        if (!allAnnots.isEmpty()){
            status.info("\tAnnotations for QTLs and Variants being made: "+allAnnots.size());
            dao.insertAnnotationsBatch(allAnnots);
        }
        if (!qtlRgdIds.isEmpty()){
            status.info("\tNew rgd_ref_rgd objects being made: " + qtlRgdIds.size());
            dao.insertRgdRefRgd(refKey,qtlRgdIds);
        }

    }


    XdbId createXdb(QTL gwasQtl) throws Exception{
        XdbId x = new XdbId();
        x.setAccId("P50DA037844");
        x.setLinkText("P50DA037844");
        x.setRgdId(gwasQtl.getRgdId());
        Date date = new Date();
        x.setCreationDate(date);
        x.setModificationDate(date);
        x.setSrcPipeline("GWAS Catalog");
        x.setXdbKey(162);
        return x;
    }

    boolean checkAnnotationExist(int annotRgdId, Term term) throws Exception{
        List<Annotation> annots = dao.getAnnotations(annotRgdId, term.getAccId(), getCreatedBy());
        if (term.isObsolete()){
            dao.deleteAnnotations(annots);
            return term.isObsolete();
        }
        if (!annots.isEmpty()){
            // update last modified data
            dao.updateLastModifiedAnnots(annots);
        }
        return !annots.isEmpty(); // if none, false
    }

    boolean checkRefAssocExist(int rgdId) throws Exception {
        List<Reference> refs = dao.getReferenceAssociations(rgdId);
        for (Reference ref : refs){
            if (ref.getKey()==refKey)
                return true;
        }
        return false;
    }

    void removeStaleAnnots()throws Exception{
        Date date0 = new Date();
        long time0 = System.currentTimeMillis();
        logStatus.info("   "+dao.getConnectionInfo());

        SimpleDateFormat sdt = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        logStatus.info("GWAS Annotation Pipeline started at "+sdt.format(date0));

        Date dtStart = Utils.addDaysToDate(new Date(), -2);
        String[] aspects = {"L","V"};
        for (String aspect : aspects) {
            deleteObsoleteAnnotations(getCreatedBy(), dtStart, getDeleteThresholdForStaleAnnotations(), getRefRgdId(), "RAT_GWAS",aspect);
        }
        logStatus.info("\nTotal pipeline runtime -- elapsed time: "+
                Utils.formatElapsedTime(time0,System.currentTimeMillis()));
    }

    int deleteObsoleteAnnotations(int createdBy, Date dt, String staleAnnotDeleteThresholdStr, int refRgdId, String dataSource, String aspect) throws Exception{

        // convert delete-threshold string to number; i.e. '5%' --> '5'
        int staleAnnotDeleteThresholdPerc = Integer.parseInt(staleAnnotDeleteThresholdStr.substring(0, staleAnnotDeleteThresholdStr.length()-1));
        // compute maximum allowed number of stale annots to be deleted
        int annotCount = dao.getCountOfAnnotationsByReference(refRgdId, dataSource, aspect);
        int staleAnnotDeleteLimit = (staleAnnotDeleteThresholdPerc * annotCount) / 100;

        List<Annotation> staleAnnots = dao.getAnnotationsModifiedBeforeTimestamp(createdBy, dt, aspect);

        logStatus.info("ANNOTATIONS_COUNT: "+annotCount);
        if( staleAnnots.size()> 0 ) {
            logStatus.info("   stale annotation delete limit (" + staleAnnotDeleteThresholdStr + "): " + staleAnnotDeleteLimit);
            logStatus.info("   stale annotations to be deleted: " + staleAnnots.size());
        }

        if( staleAnnots.size()> staleAnnotDeleteLimit ) {
            logStatus.warn("*** DELETE of stale annots aborted! *** "+staleAnnotDeleteThresholdStr+" delete threshold exceeded!");
            return 0;
        }

//        List<Integer> staleAnnotKeys = new ArrayList<>();
//        for( Annotation ann: staleAnnots ) {
////            logAnnotsDeleted.debug("DELETE "+ann.dump("|"));
//            staleAnnotKeys.add(ann.getKey());
//        }
        dao.deleteAnnotations(staleAnnots);
        return staleAnnots.size();
    }

    public void setVersion(String version) {
        this.version = version;
    }

    public String getVersion() {
        return version;
    }

    public void setCreatedBy(int createdBy) {
        this.createdBy = createdBy;
    }

    public int getCreatedBy() {
        return createdBy;
    }

    public void setRefRgdId(int refRgdId) {
        this.refRgdId = refRgdId;
    }

    public int getRefRgdId() {
        return refRgdId;
    }

    public void setRefKey(int refKey) {
        this.refKey = refKey;
    }

    public int getRefKey() {
        return refKey;
    }

    public void setDeleteThresholdForStaleAnnotations(String deleteThresholdForStaleAnnotations) {
        this.deleteThresholdForStaleAnnotations = deleteThresholdForStaleAnnotations;
    }

    public String getDeleteThresholdForStaleAnnotations() {
        return deleteThresholdForStaleAnnotations;
    }
}
