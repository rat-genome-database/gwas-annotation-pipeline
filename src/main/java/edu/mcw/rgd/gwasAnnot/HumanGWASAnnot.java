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
import java.util.Map;

public class HumanGWASAnnot {
    Logger logStatus = LogManager.getLogger("status");
    Logger logDeleteAnnots = LogManager.getLogger("deleteAnnots");
    Logger obsoleteEfo = LogManager.getLogger("obsoleteEfo");
    DAO dao = new DAO();
    String version;
    int refRgdId;
    int createdBy;
    int refKey;
    String deleteThresholdForStaleAnnotations;
    HashMap<String, Integer> ontTermNotFoundFreqMap;

    HashMap<String, String> ontIdToAspectMap;

    void createQtlAnnots() throws Exception{
//        int qtlNum = 1;
        long time0 = System.currentTimeMillis();
        Date date0 = new Date();

        logStatus.info("   "+dao.getConnectionInfo());

        SimpleDateFormat sdt = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        logStatus.info("GWAS Annotation Pipeline started at "+sdt.format(date0));

        logStatus.info("\tStarting run for GWAS QTLs");

        ontIdToAspectMap = dao.loadOntIdToAspectMap();

        ontTermNotFoundFreqMap = new HashMap<>();

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
        int obsoleteTerms = 0;
        int nullTerms = 0;
        for (GWASCatalog gc : gwas){
            if (gc.getEfoId()==null)
                continue;
            if (!gc.getSnps().startsWith("rs"))
                continue;
            if (gc.getSnps().contains(" x ") || gc.getSnps().contains(",") || gc.getSnps().contains("/"))
                continue;
            if (gc.getChr()==null || gc.getpValStr()==null)
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
            else if (gc.getChr()!=null && gc.getpVal()!=null && gc.getpValMlog()!=null && gc.getSnps()!=null){
                gwasQtl = dao.getQtlByChrPValPeakRs(gc.getChr(),gc.getpVal().toString(),gc.getpValMlog().toString(),gc.getSnps());
                if (gwasQtl != null) {
                    qtlHashMap.put(gwasQtl.getPeakRsId() + "|" + gc.getpValMlog(), gwasQtl);
                    existingQtl.add(gwasQtl);
                    gc.setQtlRgdId(gwasQtl.getRgdId());
                }
            }
            if (gwasQtl == null || gwasQtl.getRgdId() == 0) {
                gwasQtl = new QTL();
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
                String efoId;
                if (eId.startsWith("OBA_VT")) {
                    // OBA_VT ids carry a VT id directly, e.g. OBA_VT0001253 -> VT:0001253
                    efoId = "VT:" + eId.substring("OBA_VT".length());
                } else {
                    efoId = eId.replace("_", ":");
                    if (!efoId.startsWith("EFO") && !efoId.startsWith("MONDO") && !efoId.startsWith("GO") && !efoId.startsWith("HP"))
                        efoId = "EFO:" + efoId;
                }
                Term t = dao.getTermByAccId(efoId);
                if (t == null && !efoId.startsWith("MONDO")) { // figure why it is null
                    //status.info("\tOnt Term not found: "+efoId);
                    addOntTermToNotFoundMap(efoId);
                    efoId = "EFO:" + efoId;
                    t = dao.getTermByAccId(efoId);
                    if (t==null) {
                        nullTerms++;
                        continue;
                    }
                }

                String aspect = ontIdToAspectMap.get(t.getOntologyId());

                String notes;
                if (efoId.startsWith("EFO"))
                    notes = "Based on the EFO term ID";
                else if (efoId.startsWith("MONDO"))
                    notes = "Based on the MONDO term ID from GWAS";
                else if (efoId.startsWith("GO"))
                    notes = "Based on the GO term ID from GWAS";
                else if (efoId.startsWith("HP"))
                    notes = "Based on the HP term ID from GWAS";
                else if (efoId.startsWith("VT"))
                    notes = "Based on the VT term ID from GWAS";
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
                    a.setAspect(aspect);
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
                        Term term = dao.getTermByAccId(ts.getTermAcc());
                        if (term == null) {
                            nullTerms++;
                            continue;
                        }
                        if (term.isObsolete()) {
                            obsoleteTerms++;
                            continue;
                        }
                        annot.setAspect(ontIdToAspectMap.get(term.getOntologyId()));

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
            logStatus.info("\tGWAS QTLs already existing: "+existingQtl.size());
        if (!updateName.isEmpty()){
            logStatus.info("\tGWAS QTLs having their name updated: "+updateName.size());
            dao.updateQTLNameBatch(updateName);
        }
        if (!newQtls.isEmpty()) {
            logStatus.info("\tNew QTLs being made for GWAS: "+newQtls.size());
            dao.insertQTLBatch(newQtls);
            dao.updateGwasQtlRgdIdBatch(gwas);
        }
        if (!newXdbs.isEmpty())
        {
            logStatus.info("\tNew XdbIds being made for QTLs: "+newXdbs.size());
            dao.insertGwasXdbs(newXdbs);
        }
        if (!allNotes.isEmpty()){
            logStatus.info("\tNotes being made for QTLs: "+allNotes.size());
            dao.updateNote(allNotes);
        }
        if (!allAnnots.isEmpty()){
            logStatus.info("\tAnnotations for QTLs being made: "+allAnnots.size());
            dao.insertAnnotationsBatch(allAnnots);
        }
        if (!qtlRgdIds.isEmpty()){
            logStatus.info("\tNew rgd_ref_rgd objects being made: " + qtlRgdIds.size());
            dao.insertRgdRefRgd(refKey,qtlRgdIds);
        }

        logStatus.info("\tTerms that are null: "+nullTerms);
        logStatus.info("\tTerms that are obsolete: "+ obsoleteTerms);

        dumpOntTermNotFoundMap();

        logStatus.info("\tEnding run for GWAS QTLs");
        logStatus.info("\nTotal pipeline runtime -- elapsed time: "+ Utils.formatElapsedTime(time0,System.currentTimeMillis()));
    }

    void addOntTermToNotFoundMap( String termAcc ) {

        Integer count = ontTermNotFoundFreqMap.get(termAcc);
        if( count==null ) {
            count = 0;
        }
        ontTermNotFoundFreqMap.put(termAcc, count+1);
    }

    void dumpOntTermNotFoundMap() {
        if( !ontTermNotFoundFreqMap.isEmpty() ) {
            logStatus.info("\tOnt Term not found: "+ontTermNotFoundFreqMap.size());
            for( HashMap.Entry<String,Integer> entry: ontTermNotFoundFreqMap.entrySet() ) {
                logStatus.info("\t\t"+entry.getKey()+": "+entry.getValue());
            }
        }
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

        logStatus.info("   "+dao.getConnectionInfo());

        SimpleDateFormat sdt = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        logStatus.info("GWAS Annotation Pipeline started at "+sdt.format(date0));

        logStatus.info("\tStarting run for GWAS Variants");

        ontTermNotFoundFreqMap = new HashMap<>();

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
                String efoId;
                if (eid.startsWith("OBA_VT")) {
                    // OBA_VT ids carry a VT id directly, e.g. OBA_VT0001253 -> VT:0001253
                    efoId = "VT:" + eid.substring("OBA_VT".length());
                } else {
                    efoId = eid.replace("_", ":");
                    if( !efoId.startsWith("EFO") ) {
                        efoId = "EFO:" + efoId;
                    }
                }
                Term t = dao.getTermByAccId(efoId);
                if( t == null ) { // figure why it is null
                    efoId = "EFO:" + efoId;
                    t = dao.getTermByAccId(efoId);
                    if( t==null ) {
                        addOntTermToNotFoundMap(efoId);
                        continue;
                    }
                }
                String notes;
                if (t.getAccId().startsWith("EFO"))
                    notes = "Based on the EFO term ID";
                else if (t.getAccId().startsWith("MONDO"))
                    notes = "Based on the MONDO term ID from GWAS";
                else if (t.getAccId().startsWith("GO"))
                    notes = "Based on the GO term ID from GWAS";
                else if (t.getAccId().startsWith("HP"))
                    notes = "Based on the HP term ID from GWAS";
                else if (t.getAccId().startsWith("VT"))
                    notes = "Based on the VT term ID from GWAS";
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
                    a.setAspect(ontIdToAspectMap.get(t.getOntologyId()));
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
                        Term term = dao.getTermByAccId(ts.getTermAcc());
                        if (term == null)
                            continue;
                        if (term.isObsolete())
                            continue;
                        annot.setAspect(ontIdToAspectMap.get(t.getOntologyId()));
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

        dumpOntTermNotFoundMap();

        if (!allAnnots.isEmpty()){
            logStatus.info("\tAnnotations being made for Variants: "+allAnnots.size());
            dao.insertAnnotationsBatch(allAnnots);
        }
        logStatus.info("\tEnding run for GWAS Variants");

        logStatus.info("\nTotal pipeline runtime -- elapsed time: "+ Utils.formatElapsedTime(time0,System.currentTimeMillis()));
    }

    void removeStaleAnnots()throws Exception{
        long time0 = System.currentTimeMillis();
        Date date0 = new Date(time0);
        logStatus.info("   "+dao.getConnectionInfo());

        SimpleDateFormat sdt = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        String msg = "GWAS Annotation Pipeline, module -removeStaleAnnots, started at "+sdt.format(date0);
        logStatus.info(msg);
        logDeleteAnnots.info(msg);

        Date dtStart = Utils.addDaysToDate(new Date(), -14);
        int deleted = deleteObsoleteAnnotations(getCreatedBy(), dtStart, getDeleteThresholdForStaleAnnotations(), getRefRgdId(), SpeciesType.HUMAN);
        msg = "Deleted "+deleted+" annotations";
        logStatus.info(msg);
        logDeleteAnnots.info(msg);

        msg = "=== Pipeline elapsed time: " + Utils.formatElapsedTime(time0,System.currentTimeMillis())+"\n";
        logStatus.info(msg);
        logDeleteAnnots.info(msg);
    }

    void updateAnnotations() throws Exception {
        logStatus.info("\tUpdating With Info field start");
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
            logStatus.info("\t\tAnnotations being updated: " + updateWith.size());
            dao.updateAnnotations(updateWith);
        }
        logStatus.info("\tUpdating With Info field end");
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

    int deleteObsoleteAnnotations(int createdBy, Date dt, String staleAnnotDeleteThresholdStr, int refRgdId, int speciesTypeKey) throws Exception {

        // convert delete-threshold string to number; i.e. '5%' --> '5'
        int staleAnnotDeleteThresholdPerc = Integer.parseInt(staleAnnotDeleteThresholdStr.substring(0, staleAnnotDeleteThresholdStr.length()-1));
        // compute maximum allowed number of stale annots to be deleted
        int annotCount = dao.getCountOfAnnotationsByReference(refRgdId, speciesTypeKey);
        int staleAnnotDeleteLimit = (staleAnnotDeleteThresholdPerc * annotCount) / 100;

        List<Annotation> staleAnnots = dao.getAnnotationsModifiedBeforeTimestamp(createdBy, dt, refRgdId, speciesTypeKey);

        logDeleteAnnots.info("\tANNOTATIONS_COUNT: "+annotCount);
        if( staleAnnots.size()> 0 ) {
            logDeleteAnnots.info("\t\tstale annotation delete limit (" + staleAnnotDeleteThresholdStr + "): " + staleAnnotDeleteLimit);
            logDeleteAnnots.info("\t\tstale annotations to be deleted: " + staleAnnots.size());
        }

        if( staleAnnots.size()> staleAnnotDeleteLimit ) {
            String msg = "*** DELETE of stale annotations aborted! *** "+staleAnnotDeleteThresholdStr+" delete threshold exceeded!\n"
                    +"  annotation count: "+annotCount+"\n"
                    +"  stale annotation delete limit (" + staleAnnotDeleteThresholdStr + "): " + staleAnnotDeleteLimit +"\n"
                    +"  stale annotations to be deleted: " + staleAnnots.size();
            logStatus.warn(msg);
            logDeleteAnnots.warn(msg);
            return 0;
        }

        // before deleting, count annotations per aspect
        HashMap<String, Integer> deletedAnnotsPerAspect = new HashMap<>();
        for( Annotation a: staleAnnots ) {
            Integer count = deletedAnnotsPerAspect.get(a.getAspect());
            if( count == null ) {
                count = 1;
            } else {
                count++;
            }
            deletedAnnotsPerAspect.put(a.getAspect(), count);
        }

        dao.deleteAnnotations(staleAnnots);

        for( Map.Entry<String, Integer> entry: deletedAnnotsPerAspect.entrySet() ) {
            logDeleteAnnots.info("    stale annots deleted for aspect "+entry.getKey()+": "+entry.getValue());
        }

        return staleAnnots.size();
    }

    XdbId createXdb(GWASCatalog g, QTL gwasQtl) {
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
