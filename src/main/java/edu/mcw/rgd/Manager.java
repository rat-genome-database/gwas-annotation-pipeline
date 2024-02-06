package edu.mcw.rgd;

import edu.mcw.rgd.datamodel.GWASCatalog;
import edu.mcw.rgd.datamodel.QTL;
import edu.mcw.rgd.datamodel.RgdId;
import edu.mcw.rgd.datamodel.ontology.Annotation;
import edu.mcw.rgd.datamodel.ontologyx.Term;
import edu.mcw.rgd.datamodel.ontologyx.TermSynonym;
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
    int createdBy;

    public static void main(String[] args) throws Exception{

        DefaultListableBeanFactory bf = new DefaultListableBeanFactory();
        new XmlBeanDefinitionReader(bf).loadBeanDefinitions(new FileSystemResource("properties/AppConfigure.xml"));
        Manager manager = (Manager) (bf.getBean("hrdpVariants"));

        try{
            manager.run();
        }catch (Exception e){
            Utils.printStackTrace(e,manager.status);
        }
    }

    public void run() throws Exception{
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
        createQtlAnnots(gwas);
        // create QTLs based on the variants, and them make annotations for them
        status.info("Total pipeline runtime -- elapsed time: "+
                Utils.formatElapsedTime(time0,System.currentTimeMillis()));

    }

    void createQtlAnnots(List<GWASCatalog> gwas) throws Exception{
//        int qtlNum = 1;
        HashMap<String, QTL> qtlHashMap = new HashMap<>(); // rsId + PVal is key to help prevent creating duplicates
        HashMap<String, List<String>> qtlToTerm = new HashMap<>(); // make sure I do not make duplicates of Annots
        List<QTL> existingQtl = new ArrayList<>();
        List<QTL> newQtls = new ArrayList<>();
        List<Annotation> allAnnots = new ArrayList<>();
        for (GWASCatalog gc : gwas){

            if (gc.getEfoId()==null)
                continue;
            if (!gc.getSnps().startsWith("rs"))
                continue;
            if (gc.getSnps().contains(" x "))
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
                gwasQtl.setName(gc.getMapTrait() + " " + "GWAS" + qtlNum + " (human)");
                gwasQtl.setChromosome(gc.getChr());
                RgdId r = dao.createRgdId(RgdId.OBJECT_KEY_QTLS, "ACTIVE", "created by GWAS annotation Pipeline", 38);
                gwasQtl.setRgdId(r.getRgdId());
//                gwasQtl.setRgdId(0);
                gwasQtl.setPeakRsId(gc.getSnps());
                gc.setQtlRgdId(r.getRgdId());
                newQtls.add(gwasQtl);
                qtlHashMap.put(gwasQtl.getPeakRsId() + "|" + gc.getpValMlog(), gwasQtl);
            }

//            qtlNum++; // not needed with sequence

            // use EFO to find term
            // RDO, CMO, VT
            // Aspects: EFO - T
            //          CMO - L
            //          RDO - D
            //          VT  - V
            String[] efoIds = gc.getEfoId().split(", ");
            if (efoIds.length>1) // could be later, an annotation is made for both EFO ids with a loop
                continue;
            Term t = new Term();
            String efoId = gc.getEfoId().replace("_",":");
            if (!efoId.startsWith("EFO"))
                efoId = "EFO:"+efoId;
            t = dao.getTermByAccId(efoId);
            if (t==null){ // figure why it is null
                continue;
            }

            qtlToTerm.computeIfAbsent(gwasQtl.getSymbol(), k -> new ArrayList<>());
            List<String> terms = qtlToTerm.get(gwasQtl.getSymbol());

            if (!checkAnnotationExist(gwasQtl.getRgdId(),efoId) && !terms.contains(t.getAccId()) ) // does not exist
            {
                Annotation a = new Annotation();
                a.setCreatedBy(getCreatedBy());
                a.setLastModifiedBy(getCreatedBy());
                a.setAnnotatedObjectRgdId(gwasQtl.getRgdId());
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
                a.setRgdObjectKey(6); // 6 - qtls, 24 - variants
                allAnnots.add(a);
                terms.add(t.getAccId());
            }

            // get TermSynonym by EFO acc

            List<TermSynonym> synonyms = dao.getTermSynonymsBySynonymName(t.getAccId());
            if (synonyms!=null && !synonyms.isEmpty()){
                // loop through and check if they match above Aspects
                for (TermSynonym ts : synonyms){
                    Annotation annot = new Annotation();
                    if (ts.getTermAcc().startsWith("CMO")){
                        annot.setAspect("L");
                    }
                    else if (ts.getTermAcc().startsWith("DOID")){
                        annot.setAspect("D");
                    }
                    else if (ts.getTermAcc().startsWith("VT")){
                        annot.setAspect("V");
                    }
                    else
                        continue;
                    Term term = dao.getTermByAccId(ts.getTermAcc());
                    if (term==null)
                        continue;
                    if (!checkAnnotationExist(gwasQtl.getRgdId(),term.getAccId()) && !terms.contains(term.getAccId()) ) {
                        annot.setCreatedBy(getCreatedBy());
                        annot.setLastModifiedBy(getCreatedBy());
                        annot.setAnnotatedObjectRgdId(gwasQtl.getRgdId());
                        annot.setCreatedDate(new Date());
                        annot.setLastModifiedDate(annot.getLastModifiedDate());
                        annot.setWithInfo("RGD:"+gwasQtl.getRgdId());
                        annot.setTerm(term.getTerm());
                        annot.setTermAcc(term.getAccId());
                        annot.setObjectName(gwasQtl.getName());
                        annot.setObjectSymbol(gwasQtl.getSymbol());
                        annot.setSpeciesTypeKey(1);
                        annot.setDataSrc("GWAS_CATALOG");
                        annot.setEvidence("IAGP");
                        annot.setRgdObjectKey(6);
                        allAnnots.add(annot);
                        terms.add(term.getAccId());
                    }
                } // end synonym for
            }
            qtlToTerm.put(gwasQtl.getSymbol(),terms);
        }

        // insert annotations, qtls,update qwas
        if (!existingQtl.isEmpty())
            status.info("\tGWAS QTLs already existing: "+existingQtl.size());
        if (!newQtls.isEmpty()) {
            status.info("\tNew QTLs being made for GWAS: "+newQtls.size());
            dao.insertQTLBatch(newQtls);
        }
        if (!allAnnots.isEmpty()){
            status.info("\tAnnotations for QTLs being made: "+allAnnots.size());
            dao.insertAnnotationsBatch(allAnnots);
        }
        dao.updateGwasQtlRgdIdBatch(gwas);
        return;
    }

    boolean checkAnnotationExist(int annotRgdId, String accId) throws Exception{
        List<Annotation> annots = dao.getAnnotations(annotRgdId, accId, getCreatedBy());
        return !annots.isEmpty(); // if none, false
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
}
