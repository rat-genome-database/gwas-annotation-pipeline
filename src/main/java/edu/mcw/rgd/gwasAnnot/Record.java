package edu.mcw.rgd.gwasAnnot;

import edu.mcw.rgd.datamodel.GWASCatalog;
import edu.mcw.rgd.datamodel.ontology.Annotation;
import edu.mcw.rgd.datamodel.ontologyx.Term;

import java.util.ArrayList;
import java.util.List;

public class Record {
    private GWASCatalog var;

    private Term term;
    private List<Annotation> varAnnotsIncomming = new ArrayList<>();
    private List<Annotation> qtlAnnotsIncomming = new ArrayList<>();
}
