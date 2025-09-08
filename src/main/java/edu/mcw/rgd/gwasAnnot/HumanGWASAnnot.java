package edu.mcw.rgd.gwasAnnot;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.datamodel.ontology.Annotation;
import edu.mcw.rgd.datamodel.ontologyx.Term;
import edu.mcw.rgd.datamodel.ontologyx.TermSynonym;
import edu.mcw.rgd.datamodel.variants.VariantMapData;
import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.math.BigDecimal;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;

public class HumanGWASAnnot {
    Logger status = LogManager.getLogger("status");
    Logger logStatus = LogManager.getLogger("deleteAnnots");
    Logger obsoleteEfo = LogManager.getLogger("obsoleteEfo");
    DAO dao = new DAO();
    String version;
    int refRgdId;
    int createdBy;
    int refKey;
    String deleteThresholdForStaleAnnotations;
    void createQtlAnnots() throws Exception{
//        int qtlNum = 1;
        long time0 = System.currentTimeMillis();
        Date date0 = new Date();

        status.info("   "+dao.getConnectionInfo());

        SimpleDateFormat sdt = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        status.info("GWAS Annotation Pipeline started at "+sdt.format(date0));

        status.info("\tStarting run for GWAS QTLs");

        List<GWASCatalog> gwas = dao.getGWASByMapKey(38);

        HashMap<String, QTL> qtlHashMap = new HashMap<>(); // rsId + PVal is key to help prevent creating duplicates
        HashMap<String, List<String>> qtlToTerm = new HashMap<>(); // make sure I do not make duplicates of Annots
        List<QTL> existingQtl = new ArrayList<>();
        List<QTL> newQtls = new ArrayList<>();
        List<QTL> updateName = new ArrayList<>();
        List<Annotation> allAnnots = new ArrayList<>();
        List<Note> allNotes = new ArrayList<>();
        List<XdbId> newXdbs = new ArrayList<>();
        List<Integer> qtlRgdIds = new ArrayList<>();
        for (GWASCatalog gc : gwas){
            if (gc.getEfoId()==null)
                continue;
            if (!gc.getSnps().startsWith("rs"))
                continue;
            if (gc.getSnps().contains(" x ") || gc.getSnps().contains(",") || gc.getSnps().contains("/"))
                continue;
            if (gc.getChr()==null)
                continue;


            // find trait based on EFO, if none, use EFO
            // make qtl for each p Val
            QTL gwasQtl = new QTL(); // need a check if QTL exists
            if (qtlHashMap.get(gc.getSnps()+"|"+gc.getpValMlog())!=null) {
                gwasQtl = qtlHashMap.get(gc.getSnps()+"|"+gc.getpValMlog());
            }
            else if (gc.getQtlRgdId() != null  && !gc.getQtlRgdId().equals(0)){
                gwasQtl = dao.getQtlByRgdId(gc.getQtlRgdId());
                qtlHashMap.put(gwasQtl.getPeakRsId() + "|" + gc.getpValMlog(), gwasQtl);
                existingQtl.add(gwasQtl);
                String symbolWOEnd = gwasQtl.getSymbol().replace("_H","");
                if (!Utils.stringsAreEqual(gwasQtl.getName(),gc.getMapTrait() + " QTL " + symbolWOEnd + " (human)")){
                    gwasQtl.setName(gc.getMapTrait() + " QTL " + gwasQtl.getSymbol() + " (human)");
                    updateName.add(gwasQtl);
                }
            }
            else {
                BigDecimal pval = gc.getpVal();
                int scale = pval.scale();
                if (scale>83){
                    gwasQtl.setPValue(0.0);
                    gwasQtl.setpValueMlog(gc.getpValMlog());
                }
                else {
                    Double pVal = gc.getpVal().doubleValue();
                    gwasQtl.setPValue(pVal);
                }

                int qtlNum = dao.GenerateNextQTLSeqForGwas();
                gwasQtl.setSymbol("GWAS" + qtlNum + "_H");
                gwasQtl.setName(gc.getMapTrait() + " QTL " + "GWAS" + qtlNum + " (human)");
                gwasQtl.setChromosome(gc.getChr());
                RgdId r = dao.createRgdId(RgdId.OBJECT_KEY_QTLS, "ACTIVE", "created by GWAS annotation Pipeline", 38);
                gwasQtl.setRgdId(r.getRgdId());
//                gwasQtl.setRgdId(0);
                gwasQtl.setPeakRsId(gc.getSnps());
                gc.setQtlRgdId(r.getRgdId());
                newQtls.add(gwasQtl);
                qtlRgdIds.add(r.getRgdId());
                qtlHashMap.put(gwasQtl.getPeakRsId() + "|" + gc.getpValMlog(), gwasQtl);
            }
            List<XdbId> xdbs = dao.getGwasXdbs(gwasQtl.getRgdId());
            if (xdbs.isEmpty()){
                XdbId x = createXdb(gc, gwasQtl);
                newXdbs.add(x);
            }

            List<Note> noteList = dao.getQTLNoteTraits(gwasQtl.getRgdId());
            if (noteList.isEmpty()){
                Note n = new Note();
                n.setRgdId(gwasQtl.getRgdId());
                n.setPublicYN("N");
                n.setNotesTypeName("qtl_trait");
                n.setNotes(gc.getMapTrait());
                allNotes.add(n);
            }
//            qtlNum++; // not needed with sequence

            // use EFO to find term
            // RDO, CMO, VT
            // Aspects: EFO - T
            //          CMO - L
            //          RDO - D
            //          VT  - V
            String[] efoIds = gc.getEfoId().split(", ");
//            if (efoIds.length>1) // could be later, an annotation is made for both EFO ids with a loop
//                continue;
            for (String eId : efoIds) {
                Term t = new Term();
                String efoId = eId.replace("_", ":");
                if (!efoId.startsWith("EFO") && !efoId.startsWith("MONDO") && !efoId.startsWith("GO") && !efoId.startsWith("HP"))
                    efoId = "EFO:" + efoId;
                t = dao.getTermByAccId(efoId);
                if (t == null && !efoId.startsWith("MONDO")) { // figure why it is null
                    status.info("\tOnt Term not found: "+efoId);
                    efoId = "EFO:" + efoId;
                    t = dao.getTermByAccId(efoId);
                    if (t==null)
                        continue;
                }
                String notes = "";
                if (efoId.startsWith("EFO"))
                    notes = "Based on the EFO term ID";
                else if (efoId.startsWith("MONDO"))
                    notes = "Based on the MONDO term ID from GWAS";
                else if (efoId.startsWith("GO"))
                    notes = "Based on the GO term ID from GWAS";
                else if (efoId.startsWith("HP"))
                    notes = "Based on the HP term ID from GWAS";
                else
                    notes = "Based on the EFO term ID";
                qtlToTerm.computeIfAbsent(gwasQtl.getSymbol(), k -> new ArrayList<>());
                List<String> terms = qtlToTerm.get(gwasQtl.getSymbol());

                if (t != null && !checkAnnotationExist(gwasQtl.getRgdId(), t) && !terms.contains(t.getAccId())) // does not exist
                {
                    Annotation a = new Annotation();
                    a.setCreatedBy(getCreatedBy());
                    a.setLastModifiedBy(getCreatedBy());
                    a.setAnnotatedObjectRgdId(gwasQtl.getRgdId());
                    a.setRefRgdId(refRgdId);
                    a.setAspect("T");
                    a.setCreatedDate(new Date());
                    a.setLastModifiedDate(a.getCreatedDate());
                    a.setTerm(t.getTerm());
                    a.setTermAcc(t.getAccId());
                    a.setWithInfo(gwasQtl.getPeakRsId());
                    a.setObjectSymbol(gwasQtl.getSymbol());
                    a.setObjectName(gwasQtl.getName());
                    a.setSpeciesTypeKey(1);
                    a.setDataSrc("GWAS_CATALOG");
                    a.setEvidence("IAGP");
                    a.setRgdObjectKey(6); // 6 - qtls
                    a.setXrefSource(gc.getPmid());
                    allAnnots.add(a);
                    terms.add(t.getAccId());
                }

                List<TermSynonym> synonyms = dao.getTermSynonymsBySynonymName(efoId);
                if (synonyms != null && !synonyms.isEmpty()) {
                    // loop through and check if they match above Aspects
                    for (TermSynonym ts : synonyms) {
                        Annotation annot = new Annotation();
                        if (ts.getTermAcc().startsWith("CMO")) {
                            annot.setAspect("L");
                        } else if (ts.getTermAcc().startsWith("DOID")) {
                            annot.setAspect("D");
                        } else if (ts.getTermAcc().startsWith("VT")) {
                            annot.setAspect("V");
                        } else if (ts.getTermAcc().startsWith("HP")) {
                            annot.setAspect("H");
                        } else
                            continue;
                        Term term = dao.getTermByAccId(ts.getTermAcc());
                        if (term == null)
                            continue;
                        if (term.isObsolete())
                            continue;
                        if (!checkAnnotationExistWithEFO(gwasQtl.getRgdId(), term, t) && !terms.contains(term.getAccId())) {
                            annot.setCreatedBy(getCreatedBy());
                            annot.setLastModifiedBy(getCreatedBy());
                            annot.setAnnotatedObjectRgdId(gwasQtl.getRgdId());
                            annot.setCreatedDate(new Date());
                            annot.setRefRgdId(refRgdId);
                            annot.setLastModifiedDate(annot.getLastModifiedDate());
                            annot.setWithInfo(gwasQtl.getPeakRsId()); // change to rsId and then deal with the rgdweb link to work with rsId
                            annot.setTerm(term.getTerm());
                            annot.setTermAcc(term.getAccId());
                            annot.setObjectName(gwasQtl.getName());
                            annot.setObjectSymbol(gwasQtl.getSymbol());
                            annot.setSpeciesTypeKey(1);
                            annot.setDataSrc("GWAS_CATALOG");
                            annot.setEvidence("IAGP");
                            annot.setRgdObjectKey(6);
                            annot.setXrefSource(gc.getPmid());
                            annot.setNotes(notes);
                            allAnnots.add(annot);
                            terms.add(term.getAccId());
                        }
                    } // end synonym for
                }

                qtlToTerm.put(gwasQtl.getSymbol(), terms);
            }
            if (!qtlRgdIds.contains(gwasQtl.getRgdId()) && !checkRefAssocExist(gwasQtl.getRgdId())){
                qtlRgdIds.add(gwasQtl.getRgdId());
            }
        }

        // insert annotations, qtls,update qwas
        if (!existingQtl.isEmpty())
            status.info("\tGWAS QTLs already existing: "+existingQtl.size());
        if (!updateName.isEmpty()){
            status.info("\tGWAS QTLs having their name updated: "+updateName.size());
            dao.updateQTLNameBatch(updateName);
        }
        if (!newQtls.isEmpty()) {
            status.info("\tNew QTLs being made for GWAS: "+newQtls.size());
            dao.insertQTLBatch(newQtls);
            dao.updateGwasQtlRgdIdBatch(gwas);
        }
        if (!newXdbs.isEmpty())
        {
            status.info("\tNew XdbIds being made for QTLs: "+newXdbs.size());
            dao.insertGwasXdbs(newXdbs);
        }
        if (!allNotes.isEmpty()){
            status.info("\tNotes being made for QTLs: "+allNotes.size());
            dao.updateNote(allNotes);
        }
        if (!allAnnots.isEmpty()){
            status.info("\tAnnotations for QTLs being made: "+allAnnots.size());
            dao.insertAnnotationsBatch(allAnnots);
        }
        if (!qtlRgdIds.isEmpty()){
            status.info("\tNew rgd_ref_rgd objects being made: " + qtlRgdIds.size());
            dao.insertRgdRefRgd(refKey,qtlRgdIds);
        }
        status.info("\tEnding run for GWAS QTLs");

        status.info("\nTotal pipeline runtime -- elapsed time: "+
                Utils.formatElapsedTime(time0,System.currentTimeMillis()));
        return;
    }

    boolean checkRefAssocExist(int rgdId) throws Exception {
        List<Reference> refs = dao.getReferenceAssociations(rgdId);
        for (Reference ref : refs){
            if (ref.getKey()==refKey)
                return true;
        }
        return false;
    }

    void createVariantAnnots() throws Exception{
        long time0 = System.currentTimeMillis();
        Date date0 = new Date();

        status.info("   "+dao.getConnectionInfo());

        SimpleDateFormat sdt = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        status.info("GWAS Annotation Pipeline started at "+sdt.format(date0));

        status.info("\tStarting run for GWAS Variants");

        List<GWASCatalog> gwas = dao.getGWASByMapKey(38);

        HashMap<Long,List<String>> varToTerm = new HashMap<>();
        List<Annotation> allAnnots = new ArrayList<>();
        for (GWASCatalog gc : gwas) {
            if (gc.getEfoId() == null)
                continue;
            if (!gc.getSnps().startsWith("rs"))
                continue;
            if (gc.getSnps().contains(" x ") || gc.getSnps().contains(",") || gc.getSnps().contains("/"))
                continue;
            if (gc.getChr() == null)
                continue;

            List<VariantMapData> vars = dao.getAllActiveVariantsWithRsId(gc.getSnps());
            VariantMapData vmd = new VariantMapData();
            boolean found = false;
            for (VariantMapData v : vars){
                if (Utils.stringsAreEqual(gc.getStrongSnpRiskallele(),v.getVariantNucleotide())){
                    vmd = v;
                    found = true;
                }
            }
            if (!found)
                continue;

            int rgdId = Math.toIntExact(vmd.getId());
            String[] efoIds = gc.getEfoId().split(", ");

            for (String eid : efoIds) {
                Term t = new Term();
                String efoId =  eid.replace("_", ":");
                if (!efoId.startsWith("EFO") && !efoId.startsWith("MONDO") && !efoId.startsWith("GO") && !efoId.startsWith("HP"))
                    efoId = "EFO:" + efoId;
                t = dao.getTermByAccId(efoId);
                if (t == null) { // figure why it is null
                    status.info("\tOnt Term not found: "+efoId);
                    efoId = "EFO:" + efoId;
                    t = dao.getTermByAccId(efoId);
                    if (t==null)
                        continue;
                }
                String notes = "";
                if (t.getAccId().startsWith("EFO"))
                    notes = "Based on the EFO term ID";
                else if (t.getAccId().startsWith("MONDO"))
                    notes = "Based on the MONDO term ID from GWAS";
                else if (t.getAccId().startsWith("GO"))
                    notes = "Based on the GO term ID from GWAS";
                else if (t.getAccId().startsWith("HP"))
                    notes = "Based on the HP term ID from GWAS";
                else
                    notes = "Based on the EFO term ID";

                varToTerm.computeIfAbsent(vmd.getId(), k -> new ArrayList<>());
                List<String> terms = varToTerm.get(vmd.getId());

                if (!checkAnnotationExist(rgdId, t) && !terms.contains(t.getAccId())) {
                    Annotation a = new Annotation();
                    a.setCreatedBy(getCreatedBy());
                    a.setLastModifiedBy(getCreatedBy());
                    a.setAnnotatedObjectRgdId(rgdId);
                    a.setRefRgdId(refRgdId);
                    a.setAspect("T");
                    a.setCreatedDate(new Date());
                    a.setLastModifiedDate(a.getLastModifiedDate());
//                    a.setWithInfo("RGD:"+vmd.getId());
                    a.setTerm(t.getTerm());
                    a.setTermAcc(t.getAccId());
//                    a.setObjectName(vmd.getRsId());
                    a.setObjectSymbol(vmd.getRsId());
                    a.setSpeciesTypeKey(1);
                    a.setDataSrc("GWAS_CATALOG");
                    a.setEvidence("IAGP");
                    a.setRgdObjectKey(7); // 7 - variants
                    a.setXrefSource(gc.getPmid());
                    allAnnots.add(a);
                    terms.add(t.getAccId());
                }

                List<TermSynonym> synonyms = dao.getTermSynonymsBySynonymName(t.getAccId());
                if (synonyms != null && !synonyms.isEmpty()) {
                    // loop through and check if they match above Aspects
                    for (TermSynonym ts : synonyms) {
                        Annotation annot = new Annotation();
                        if (ts.getTermAcc().startsWith("CMO")) {
                            annot.setAspect("L");
                        } else if (ts.getTermAcc().startsWith("DOID")) {
                            annot.setAspect("D");
                        } else if (ts.getTermAcc().startsWith("VT")) {
                            annot.setAspect("V");
                        } else if(ts.getTermAcc().startsWith("HP")){
                            annot.setAspect("H");
                        } else
                            continue;
                        Term term = dao.getTermByAccId(ts.getTermAcc());
                        if (term == null)
                            continue;
                        if (term.isObsolete())
                            continue;
                        if (!checkAnnotationExistWithEFO(rgdId, term, t) && !terms.contains(term.getAccId())) {
                            annot.setCreatedBy(getCreatedBy());
                            annot.setLastModifiedBy(getCreatedBy());
                            annot.setAnnotatedObjectRgdId(rgdId);
                            annot.setCreatedDate(new Date());
                            annot.setRefRgdId(refRgdId);
                            annot.setLastModifiedDate(annot.getLastModifiedDate());
//                            annot.setWithInfo("RGD:" + vmd.getId());
                            annot.setTerm(term.getTerm());
                            annot.setTermAcc(term.getAccId());
//                            annot.setObjectName(vmd.getRsId());
                            annot.setObjectSymbol(vmd.getRsId());
                            annot.setSpeciesTypeKey(1);
                            annot.setDataSrc("GWAS_CATALOG");
                            annot.setEvidence("IAGP");
                            annot.setRgdObjectKey(7);
                            annot.setXrefSource(gc.getPmid());
                            annot.setNotes(notes);
                            allAnnots.add(annot);
                            terms.add(term.getAccId());
                        }
                    } // end synonym for
                }
                varToTerm.put(vmd.getId(),terms);
            } // end efo for
        }
        if (!allAnnots.isEmpty()){
            status.info("\tAnnotations being made for Variants: "+allAnnots.size());
            dao.insertAnnotationsBatch(allAnnots);
        }
        status.info("\tEnding run for GWAS Variants");

        status.info("\nTotal pipeline runtime -- elapsed time: "+
                Utils.formatElapsedTime(time0,System.currentTimeMillis()));
    }

    void removeStaleAnnots()throws Exception{
        Date date0 = new Date();
        long time0 = System.currentTimeMillis();
        logStatus.info("   "+dao.getConnectionInfo());

        SimpleDateFormat sdt = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        logStatus.info("GWAS Annotation Pipeline started at "+sdt.format(date0));

        Date dtStart = Utils.addDaysToDate(new Date(), -7);
        String[] aspects = {"T","L","D","V","H"};
        for (String aspect : aspects) {
            String ont = "";
            switch (aspect){
                case "T":
                    ont="EFO";
                    break;
                case "L":
                    ont="CMO";
                    break;
                case "D":
                    ont="DOID";
                    break;
                case "V":
                    ont="VT";
                    break;
                case "H":
                    ont="HP";
                    break;
            }
            logStatus.info("Running for Ontology: "+ont);
            deleteObsoleteAnnotations(getCreatedBy(), dtStart, getDeleteThresholdForStaleAnnotations(), getRefRgdId(), "GWAS_CATALOG",aspect);
        }
        logStatus.info("\nTotal pipeline runtime -- elapsed time: "+
                Utils.formatElapsedTime(time0,System.currentTimeMillis()));
    }

    void updateAnnotations() throws Exception {
        status.info("\tUpdating With Info field start");
        List<GWASCatalog> gwas = dao.getGWASByMapKey(38);
        List<Annotation> updateWith = new ArrayList<>();
        List<Integer> rgdIds = new ArrayList<>();
        for (GWASCatalog g : gwas){
            if (g.getQtlRgdId()==0)
                continue;
            if (rgdIds.contains(g.getQtlRgdId()))
                continue;
            QTL gwasQtl = dao.getQtlByRgdId(g.getQtlRgdId());
            List<Annotation> annots = dao.getAnnotations(gwasQtl.getRgdId());
            for (Annotation a : annots){
                if (a.getWithInfo()==null || a.getWithInfo().startsWith("RGD:")) {
                    a.setWithInfo(gwasQtl.getPeakRsId());
                    a.setLastModifiedBy(getCreatedBy());
                    updateWith.add(a);
                }
            }
            rgdIds.add(gwasQtl.getRgdId());
        }
        if (!updateWith.isEmpty()){
            status.info("\t\tAnnotations being updated: " + updateWith.size());
            dao.updateAnnotations(updateWith);
        }
        status.info("\tUpdating With Info field end");
    }

    void checkDbSnp()throws Exception{
        List<GWASCatalog> gwas = dao.getGWASCatalog();
        Logger snpSum = LogManager.getLogger("snpSummary");
        snpSum.info("GWAS SNPs not found in DB_SNP:");
        for (GWASCatalog gc : gwas){
            if (!gc.getSnps().startsWith("rs"))
                continue;
            if (gc.getSnps().contains(" x ") || gc.getSnps().contains(",") || gc.getSnps().contains("/"))
                continue;
            List<dbSnp> snps = dao.getSnpBySnpName(gc.getSnps());
            if (snps.isEmpty()) {
//                if (snps.size()>1)
//                    System.out.println(gc.getSnps());
////                for (dbSnp snp : snps) {
////                    // cross check allele if more than one
////
////                }
//            }
//            else {
                snpSum.info(gc.getSnps());
            }
        }
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

    boolean checkAnnotationExistWithEFO(int annotRgdId, Term term, Term efo) throws Exception{
        List<Annotation> annots = dao.getAnnotations(annotRgdId, term.getAccId(), getCreatedBy());
        if (term.isObsolete()){
            dao.deleteAnnotations(annots);
            return term.isObsolete();
        }
        if (efo != null && efo.isObsolete() && !annots.isEmpty()){
            for (Annotation a : annots) {
                obsoleteEfo.info("Annotation Based on Obsolete EFO:\n" + a.dump("|"));
            }
        }
        if (!annots.isEmpty()){
            dao.updateLastModifiedAnnots(annots);
            return true;
        }
        return false; // if none, false
    }

    int deleteObsoleteAnnotations(int createdBy, Date dt, String staleAnnotDeleteThresholdStr, int refRgdId, String dataSource, String aspect) throws Exception{

        // convert delete-threshold string to number; i.e. '5%' --> '5'
        int staleAnnotDeleteThresholdPerc = Integer.parseInt(staleAnnotDeleteThresholdStr.substring(0, staleAnnotDeleteThresholdStr.length()-1));
        // compute maximum allowed number of stale annots to be deleted
        int annotCount = dao.getCountOfAnnotationsByReference(refRgdId, dataSource, aspect);
        int staleAnnotDeleteLimit = (staleAnnotDeleteThresholdPerc * annotCount) / 100;

        List<Annotation> staleAnnots = dao.getAnnotationsModifiedBeforeTimestamp(createdBy, dt, aspect);

        logStatus.info("\tANNOTATIONS_COUNT: "+annotCount);
        if( staleAnnots.size()> 0 ) {
            logStatus.info("\t\tstale annotation delete limit (" + staleAnnotDeleteThresholdStr + "): " + staleAnnotDeleteLimit);
            logStatus.info("\t\tstale annotations to be deleted: " + staleAnnots.size());
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
//        dao.deleteAnnotations(staleAnnots);
        return staleAnnots.size();
    }

    XdbId createXdb(GWASCatalog g, QTL gwasQtl) throws Exception{
        XdbId x = new XdbId();
        x.setAccId(g.getStudyAcc());
        x.setLinkText(g.getStudyAcc());
        x.setRgdId(gwasQtl.getRgdId());
        Date date = new Date();
        x.setCreationDate(date);
        x.setModificationDate(date);
        x.setSrcPipeline("GWAS Catalog");
        x.setXdbKey(dao.getXdbKey());
        return x;
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
    public int getRefKey(){
        return refKey;
    }

    public void setDeleteThresholdForStaleAnnotations(String deleteThresholdForStaleAnnotations) {
        this.deleteThresholdForStaleAnnotations = deleteThresholdForStaleAnnotations;
    }

    public String getDeleteThresholdForStaleAnnotations() {
        return deleteThresholdForStaleAnnotations;
    }
}
