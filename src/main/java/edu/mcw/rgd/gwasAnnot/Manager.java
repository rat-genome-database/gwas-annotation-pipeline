package edu.mcw.rgd.gwasAnnot;

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

    public static void main(String[] args) throws Exception{

        DefaultListableBeanFactory bf = new DefaultListableBeanFactory();
        new XmlBeanDefinitionReader(bf).loadBeanDefinitions(new FileSystemResource("properties/AppConfigure.xml"));
        HumanGWASAnnot ha = (HumanGWASAnnot) (bf.getBean("gwasAnnot"));
        try{
//            manager.run(args);
            for (String arg : args) {
                switch (arg) {
                    case "-checkDbDnp":
                        ha.checkDbSnp();
                        break;
                    case "-qtlAnnotRun":
                        ha.createQtlAnnots();
                        break;
                    case "-varAnnotRun":
                        ha.createVariantAnnots();
                        break;
                    case "-updateAnnots":
                        ha.updateAnnotations();
                        break;
                    case "-removeStaleAnnots":
                        ha.removeStaleAnnots();
                        break;
                    case "-ratAnnotRun":
                        RatAnnot ra = (RatAnnot) (bf.getBean("ratAnnot"));
                        ra.run();
                        break;
                    case "-removeStaleRatAnnots":
                        RatAnnot ra1 = (RatAnnot) (bf.getBean("ratAnnot"));
                        ra1.removeStaleAnnots();
                        break;
                }
            }
        }catch (Exception e){
            Utils.printStackTrace(e,ha.status);
        }
    }

}
