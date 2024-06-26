package edu.mcw.rgd;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.datamodel.ontology.Annotation;
import edu.mcw.rgd.datamodel.ontologyx.Term;
import edu.mcw.rgd.datamodel.ontologyx.TermSynonym;
import edu.mcw.rgd.datamodel.variants.VariantMapData;
import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.beans.factory.support.DefaultListableBeanFactory;
import org.springframework.beans.factory.xml.XmlBeanDefinitionReader;
import org.springframework.core.io.FileSystemResource;

import java.math.BigDecimal;
import java.text.SimpleDateFormat;
import java.util.*;

public class Manager {
    Logger status = LogManager.getLogger("status");
    DAO dao = new DAO();
    String version;
    int refRgdId;
    int createdBy;

    public static void main(String[] args) throws Exception{

        DefaultListableBeanFactory bf = new DefaultListableBeanFactory();
        new XmlBeanDefinitionReader(bf).loadBeanDefinitions(new FileSystemResource("properties/AppConfigure.xml"));
        Manager manager = (Manager) (bf.getBean("gwasAnnot"));

        try{
            manager.run(args);
        }catch (Exception e){
            Utils.printStackTrace(e,manager.status);
        }
    }

    public void run(String[] args) throws Exception{
        long time0 = System.currentTimeMillis();
        Date date0 = new Date();

        status.info("   "+dao.getConnectionInfo());

        SimpleDateFormat sdt = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        status.info("GWAS Annotation Pipeline started at "+sdt.format(date0));

//        int originalAnnotCount = dao.getAnnotationsModifiedBeforeTimestamp(date0, getCreatedBy()).size();
//        if (originalAnnotCount!=0){
//            status.info("ANNOT COUNT ORIGINAL: " + Utils.formatThousands(originalAnnotCount));
//        }
        // create annotations for the variants
        List<GWASCatalog> gwas = dao.getGWASCatalog();
        for (String arg : args) {
            switch (arg) {
                case "-checkDbDnp" -> checkDbSnp();
                case "-qtlAnnotRun" -> createQtlAnnots(gwas);
                case "-varAnnotRun" -> createVariantAnnots(gwas);
            }
        }

        // create QTLs based on the variants, and them make annotations for them
        status.info("\nTotal pipeline runtime -- elapsed time: "+
                Utils.formatElapsedTime(time0,System.currentTimeMillis()));

    }

    void createQtlAnnots(List<GWASCatalog> gwas) throws Exception{
//        int qtlNum = 1;
        HashMap<String, QTL> qtlHashMap = new HashMap<>(); // rsId + PVal is key to help prevent creating duplicates
        HashMap<String, List<String>> qtlToTerm = new HashMap<>(); // make sure I do not make duplicates of Annots
        List<QTL> existingQtl = new ArrayList<>();
        List<QTL> newQtls = new ArrayList<>();
        List<QTL> updateName = new ArrayList<>();
        List<Annotation> allAnnots = new ArrayList<>();
        List<Note> allNotes = new ArrayList<>();
        List<XdbId> newXdbs = new ArrayList<>();
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
                if (!efoId.startsWith("EFO"))
                    efoId = "EFO:" + efoId;
                t = dao.getTermByAccId(efoId);
                if (t == null) { // figure why it is null
                    continue;
                }

                qtlToTerm.computeIfAbsent(gwasQtl.getSymbol(), k -> new ArrayList<>());
                List<String> terms = qtlToTerm.get(gwasQtl.getSymbol());

                if (!checkAnnotationExist(gwasQtl.getRgdId(), efoId) && !terms.contains(t.getAccId())) // does not exist
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

                // get TermSynonym by EFO acc

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
                        } else if(ts.getTermAcc().startsWith("MP")){
                            annot.setAspect("N");
                        }else if(ts.getTermAcc().startsWith("HP")){
                            annot.setAspect("H");
                        } else
                            continue;
                        Term term = dao.getTermByAccId(ts.getTermAcc());
                        if (term == null)
                            continue;
                        if (!checkAnnotationExist(gwasQtl.getRgdId(), term.getAccId()) && !terms.contains(term.getAccId())) {
                            annot.setCreatedBy(getCreatedBy());
                            annot.setLastModifiedBy(getCreatedBy());
                            annot.setAnnotatedObjectRgdId(gwasQtl.getRgdId());
                            annot.setCreatedDate(new Date());
                            annot.setRefRgdId(refRgdId);
                            annot.setLastModifiedDate(annot.getLastModifiedDate());
                            annot.setWithInfo("RGD:" + gwasQtl.getRgdId());
                            annot.setTerm(term.getTerm());
                            annot.setTermAcc(term.getAccId());
                            annot.setObjectName(gwasQtl.getName());
                            annot.setObjectSymbol(gwasQtl.getSymbol());
                            annot.setSpeciesTypeKey(1);
                            annot.setDataSrc("GWAS_CATALOG");
                            annot.setEvidence("IAGP");
                            annot.setRgdObjectKey(6);
                            annot.setXrefSource(gc.getPmid());
                            allAnnots.add(annot);
                            terms.add(term.getAccId());
                        }
                    } // end synonym for
                }
                qtlToTerm.put(gwasQtl.getSymbol(), terms);
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

        return;
    }

    void createVariantAnnots(List<GWASCatalog> catalog) throws Exception{
        HashMap<Long,List<String>> varToTerm = new HashMap<>();
        List<Annotation> allAnnots = new ArrayList<>();
        for (GWASCatalog gc : catalog) {
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
                if (!efoId.startsWith("EFO"))
                    efoId = "EFO:" + efoId;
                t = dao.getTermByAccId(efoId);
                if (t == null) { // figure why it is null
                    continue;
                }

                varToTerm.computeIfAbsent(vmd.getId(), k -> new ArrayList<>());
                List<String> terms = varToTerm.get(vmd.getId());

                if (!checkAnnotationExist(rgdId, t.getAccId()) && !terms.contains(t.getAccId())) {
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
                        } else if(ts.getTermAcc().startsWith("MP")){
                            annot.setAspect("N");
                        }else if(ts.getTermAcc().startsWith("HP")){
                            annot.setAspect("H");
                        } else
                            continue;
                        Term term = dao.getTermByAccId(ts.getTermAcc());
                        if (term == null)
                            continue;

                        if (!checkAnnotationExist(rgdId, term.getAccId()) && !terms.contains(term.getAccId())) {
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

    boolean checkAnnotationExist(int annotRgdId, String accId) throws Exception{
        List<Annotation> annots = dao.getAnnotations(annotRgdId, accId, getCreatedBy());
        return !annots.isEmpty(); // if none, false
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
}
