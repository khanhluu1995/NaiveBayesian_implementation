import java.util.Hashtable;
import java.util.Scanner;

public class NaiveBayesian_Implementation {
//    private final double nucleus = 0.384;
//    private final double cytoplasm = 0.21421997238840312;
//    private final double cytoskeleton = 0.08950759318913944;
//    private final double mitochondria = 0.0796134376438104;
//    private final double plasma_membrane = 0.07547169811320754;
//    private final double ER = 0.05499309710078233;
//    private final double golgi = 0.040266912103083294;
//    private final double vacuole = 0.0310630464795214;
//    private final double peroxisome = 0.011965025310630465 ;
//    private final double  endosome = 0.00529222273354809;
//    private final double extracellular = 0.0025310630464795212;
//    private final double integral_membrane = 9.203865623561896E-4;
//    private final double  cell_wall = 6.902899217671422E-4;
//    private final double  lipid_particles= 4.601932811780948E-4;

    // the 13 possible locals
    private  double p_nucleus  ;
    private  double p_cytoplasm  ;
    private  double p_cytoskeleton ;
    private  double p_mitochondria ;
    private  double p_plasma_membrane ;
    private  double p_ER ;
    private  double p_golgi ;
    private  double p_vacuole ;
    private  double p_peroxisome ;
    private  double p_endosome;
    private  double p_extracellular;
    private  double p_integral_membrane;
    private  double p_cell_wall ;
    private  double p_lipid_particles;
    private double p_transport_vesicles;

    private  double nucleus  ;
    private  double cytoplasm  ;
    private  double cytoskeleton ;
    private  double mitochondria ;
    private  double plasma_membrane ;
    private  double ER ;
    private  double golgi ;
    private  double vacuole ;
    private  double peroxisome ;
    private  double endosome;
    private  double extracellular;
    private  double integral_membrane;
    private  double cell_wall ;
    private  double lipid_particles;
    private double transport_vesicles;

    //the 7 attributes for calculation
    private double[] countGeneID = new double[15];
    private double[] countEssential = new double[15];
    private double[] countClass = new double[15];
    private double[] countComplex = new double[15];
    private double[]  countPhenotype = new double[15];
    private double[] countMotif = new double[15];
    private double[] countChromosome = new double[15];

    private double[] probabilityResult = new double[15];
    String[][] myDataSet;
    Hashtable<String, String> uniqueLocalization = new Hashtable<>();
    Scanner scanner;


    private double totalData;

    public NaiveBayesian_Implementation(Hashtable<String, Integer> localizationHashtable, String[][] myDataSet, double totalData) {
        this.p_nucleus = localizationHashtable.get("nucleus") / totalData;
        this.nucleus = localizationHashtable.get("nucleus");

        this.p_cytoplasm = localizationHashtable.get("cytoplasm") / totalData;
        this.cytoplasm = localizationHashtable.get("cytoplasm");

        this.p_cytoskeleton = localizationHashtable.get("cytoskeleton") / totalData;
        this.cytoskeleton = localizationHashtable.get("cytoskeleton");

        this.p_mitochondria = localizationHashtable.get("mitochondria") / totalData;
        this.mitochondria = localizationHashtable.get("mitochondria");

        this.p_plasma_membrane = localizationHashtable.get("plasma membrane") / totalData;
        this.plasma_membrane = localizationHashtable.get("plasma membrane");

        this.p_ER = localizationHashtable.get("ER")/ totalData;
        this.ER = localizationHashtable.get("ER");

        this.p_golgi = localizationHashtable.get("golgi") / totalData;
        this.golgi = localizationHashtable.get("golgi");

        this.p_vacuole = localizationHashtable.get("vacuole") / totalData;
        this.vacuole = localizationHashtable.get("vacuole");

        this.p_peroxisome = localizationHashtable.get("peroxisome") / totalData;
        this.peroxisome = localizationHashtable.get("peroxisome");

        this.p_endosome = localizationHashtable.get("endosome") / totalData;
        this.endosome = localizationHashtable.get("endosome");

        this.p_extracellular = localizationHashtable.get("extracellular") / totalData;
        this.extracellular = localizationHashtable.get("extracellular");

        this.p_integral_membrane = localizationHashtable.get("integral membrane") / totalData;
        this.integral_membrane = localizationHashtable.get("integral membrane");

        this.p_cell_wall = localizationHashtable.get("cell wall") / totalData;
        this.cell_wall = localizationHashtable.get("cell wall");

        this.p_lipid_particles = localizationHashtable.get("lipid particles") / totalData;
        this.lipid_particles = localizationHashtable.get("lipid particles");

        this.p_transport_vesicles = localizationHashtable.get("transport vesicles") /totalData;
        this.transport_vesicles = localizationHashtable.get("transport vesicles");

        this.myDataSet = myDataSet;
        this.totalData = totalData;
        setupFunction(this.countGeneID);
        setupFunction(this.countChromosome);
        setupFunction(this.countClass);
        setupFunction(this.countComplex);
        setupFunction(this.countEssential);
        setupFunction(this.countMotif);
        setupFunction(this.countPhenotype);
        setupFunction(this.probabilityResult);
    }

    private void setupFunction(double[] getSetup){
        for (int i = 0; i < getSetup.length; i++){
            getSetup[i] = 0;
        }
    }

    public double getP_nucleus() {
        return p_nucleus;
    }

    public double getP_cytoplasm() {
        return p_cytoplasm;
    }

    public double getP_cytoskeleton() {
        return p_cytoskeleton;
    }

    public double getP_mitochondria() {
        return p_mitochondria;
    }

    public double getP_plasma_membrane() {
        return p_plasma_membrane;
    }

    public double getP_ER() {
        return p_ER;
    }

    public double getP_golgi() {
        return p_golgi;
    }

    public double getP_vacuole() {
        return p_vacuole;
    }

    public double getP_peroxisome() {
        return p_peroxisome;
    }

    public double getP_endosome() {
        return p_endosome;
    }

    public double getP_extracellular() {
        return p_extracellular;
    }

    public double getP_integral_membrane() {
        return p_integral_membrane;
    }

    public double getP_cell_wall() {
        return p_cell_wall;
    }

    public double getP_lipid_particles() {
        return p_lipid_particles;
    }

    public double getTotalData() {
        return totalData;
    }

    public String predictLocalization(String[] oneLineTestData){

        findProbEssential(oneLineTestData[1]);
        findProbClass(oneLineTestData[2]);
        findProbComplex(oneLineTestData[3]);
        findProbPhenotype(oneLineTestData[4]);
        findProbMotif(oneLineTestData[5]);
        findProbChromosome(oneLineTestData[6]);






        int chosen = calculateFinalProb(probabilityResult);


        switch (chosen){
            case 1:
                return  "nucleus";


            case 2:
                return  "cytoplasm";

            case 3:
                return  "cytoskeleton";

            case 4:
                return  "mitochondria";

            case 5:
                return "plasma membrane";

            case 6:
                return "ER";

            case 7:
                return "golgi";

            case 8:
                return "vacuole";

            case 9:return "peroxisome";
            case 10: return "endosome";
            case 11: return "extracellular";
            case 12: return "integral membrane";
            case 13: return "cell wall";
            case 14: return "lipid particles";
            case 15: return "transport vesicles";

            default: return "cytoplasm";
        }

    }

    private void findProbEssential(String str){

        setupFunction(countEssential);
        for(int i = 0; i < totalData; i++){
            if(str.equals(myDataSet[i][1])  && myDataSet[i][7].equals("nucleus")){
                countEssential[0]++;
            }
            else if(str.equals(myDataSet[i][1])  && myDataSet[i][7].equals("cytoplasm")){
                countEssential[1]++;
            }
            else if(str.equals(myDataSet[i][1])  && myDataSet[i][7].equals("cytoskeleton")){
                countEssential[2]++;
            }
            else if(str.equals(myDataSet[i][1])  && myDataSet[i][7].equals("mitochondria")){
                countEssential[3]++;
            }
            else if(str.equals(myDataSet[i][1])  && myDataSet[i][7].equals("plasma membrane")){
                countEssential[4]++;
            }
            else if(str.equals(myDataSet[i][1])  && myDataSet[i][7].equals("ER")){
                countEssential[5]++;
            }
            else if(str.equals(myDataSet[i][1])  && myDataSet[i][7].equals("golgi")){
                countEssential[6]++;
            }
            else if(str.equals(myDataSet[i][1])  && myDataSet[i][7].equals("vacuole")){
                countEssential[7]++;
            }
            else if(str.equals(myDataSet[i][1])  && myDataSet[i][7].equals("peroxisome")){
                countEssential[8]++;
            }
            else if(str.equals(myDataSet[i][1])  && myDataSet[i][7].equals("endosome")){
                countEssential[9]++;
            }
            else if(str.equals(myDataSet[i][1])  && myDataSet[i][7].equals("extracellular")){
                countEssential[10]++;
            }
            else if(str.equals(myDataSet[i][1])  && myDataSet[i][7].equals("integral membrane")){
                countEssential[11]++;
            }
            else if(str.equals(myDataSet[i][1])  && myDataSet[i][7].equals("cell wall")){
                countEssential[12]++;
            }
            else if(str.equals(myDataSet[i][1])  && myDataSet[i][7].equals("lipid particles")){
                countEssential[13]++;
            }
            else if(str.equals(myDataSet[i][1])  && myDataSet[i][7].equals("transport vesicles")){
                countEssential[14]++;
            }

//            else {
//                System.out.println("input string: " + str);
//                System.out.println("data to compare string: " + myDataSet[i][1]);
//                System.out.println("name of the local at the line: " + myDataSet[i][7]);
//            }
        }
    }

    private void findProbClass(String str){

        setupFunction(countClass);
        for(int i = 0; i < totalData; i++){
            if(str.equals(myDataSet[i][2])  && myDataSet[i][7].equals("nucleus")){
                countClass[0]++;
            }
            else if(str.equals(myDataSet[i][2])  && myDataSet[i][7].equals("cytoplasm")){
                countClass[1]++;
            }
            else if(str.equals(myDataSet[i][2])  && myDataSet[i][7].equals("cytoskeleton")){
                countClass[2]++;
            }
            else if(str.equals(myDataSet[i][2])  && myDataSet[i][7].equals("mitochondria")){
                countClass[3]++;
            }
            else if(str.equals(myDataSet[i][2])  && myDataSet[i][7].equals("plasma membrane")){
                countClass[4]++;
            }
            else if(str.equals(myDataSet[i][2])  && myDataSet[i][7].equals("ER")){
                countClass[5]++;
            }
            else if(str.equals(myDataSet[i][2])  && myDataSet[i][7].equals("golgi")){
                countClass[6]++;
            }
            else if(str.equals(myDataSet[i][2])  && myDataSet[i][7].equals("vacuole")){
                countClass[7]++;
            }
            else if(str.equals(myDataSet[i][2])  && myDataSet[i][7].equals("peroxisome")){
                countClass[8]++;
            }
            else if(str.equals(myDataSet[i][2])  && myDataSet[i][7].equals("endosome")){
                countClass[9]++;
            }
            else if(str.equals(myDataSet[i][2])  && myDataSet[i][7].equals("extracellular")){
                countClass[10]++;
            }
            else if(str.equals(myDataSet[i][2])  && myDataSet[i][7].equals("integral membrane")){
                countClass[11]++;
            }
            else if(str.equals(myDataSet[i][2])  && myDataSet[i][7].equals("cell wall")){
                countClass[12]++;
            }
            else if(str.equals(myDataSet[i][2])  && myDataSet[i][7].equals("lipid particles")){
                countClass[13]++;
            }
            else if(str.equals(myDataSet[i][2])  && myDataSet[i][7].equals("transport vesicles")){
                countClass[14]++;

            }
        }
    }

    private void findProbComplex(String str){

        setupFunction(countComplex);
        for(int i = 0; i < totalData; i++){
            if(str.equals(myDataSet[i][3])  && myDataSet[i][7].equals("nucleus")){
                countComplex[0]++;
            }
            else if(str.equals(myDataSet[i][3])  && myDataSet[i][7].equals("cytoplasm")){
                countComplex[1]++;
            }
            else if(str.equals(myDataSet[i][3])  && myDataSet[i][7].equals("cytoskeleton")){
                countComplex[2]++;
            }
            else if(str.equals(myDataSet[i][3])  && myDataSet[i][7].equals("mitochondria")){
                countComplex[3]++;
            }
            else if(str.equals(myDataSet[i][3])  && myDataSet[i][7].equals("plasma membrane")){
                countComplex[4]++;
            }
            else if(str.equals(myDataSet[i][3])  && myDataSet[i][7].equals("ER")){
                countComplex[5]++;
            }
            else if(str.equals(myDataSet[i][3])  && myDataSet[i][7].equals("golgi")){
                countComplex[6]++;
            }
            else if(str.equals(myDataSet[i][3])  && myDataSet[i][7].equals("vacuole")){
                countComplex[7]++;
            }
            else if(str.equals(myDataSet[i][3])  && myDataSet[i][7].equals("peroxisome")){
                countComplex[8]++;
            }
            else if(str.equals(myDataSet[i][3])  && myDataSet[i][7].equals("endosome")){
                countComplex[9]++;
            }
            else if(str.equals(myDataSet[i][3])  && myDataSet[i][7].equals("extracellular")){
                countComplex[10]++;
            }
            else if(str.equals(myDataSet[i][3])  && myDataSet[i][7].equals("integral membrane")){
                countComplex[11]++;
            }
            else if(str.equals(myDataSet[i][3])  && myDataSet[i][7].equals("cell wall")){
                countComplex[12]++;
            }
            else if(str.equals(myDataSet[i][3])  && myDataSet[i][7].equals("lipid particles")){
                countComplex[13]++;
            }
            else if(str.equals(myDataSet[i][3])  && myDataSet[i][7].equals("transport vesicles")){
                countComplex[14]++;
            }
        }
    }

    private void findProbPhenotype(String str){

        setupFunction(countPhenotype);
        for(int i = 0; i < totalData; i++){
            if(str.equals(myDataSet[i][4])  && myDataSet[i][7].equals("nucleus")){
                countPhenotype[0]++;
            }
            else if(str.equals(myDataSet[i][4])  && myDataSet[i][7].equals("cytoplasm")){
                countPhenotype[1]++;
            }
            else if(str.equals(myDataSet[i][4])  && myDataSet[i][7].equals("cytoskeleton")){
                countPhenotype[2]++;
            }
            else if(str.equals(myDataSet[i][4])  && myDataSet[i][7].equals("mitochondria")){
                countPhenotype[3]++;
            }
            else if(str.equals(myDataSet[i][4])  && myDataSet[i][7].equals("plasma membrane")){
                countPhenotype[4]++;
            }
            else if(str.equals(myDataSet[i][4])  && myDataSet[i][7].equals("ER")){
                countPhenotype[5]++;
            }
            else if(str.equals(myDataSet[i][4])  && myDataSet[i][7].equals("golgi")){
                countPhenotype[6]++;
            }
            else if(str.equals(myDataSet[i][4])  && myDataSet[i][7].equals("vacuole")){
                countPhenotype[7]++;
            }
            else if(str.equals(myDataSet[i][4])  && myDataSet[i][7].equals("peroxisome")){
                countPhenotype[8]++;
            }
            else if(str.equals(myDataSet[i][4])  && myDataSet[i][7].equals("endosome")){
                countPhenotype[9]++;
            }
            else if(str.equals(myDataSet[i][4])  && myDataSet[i][7].equals("extracellular")){
                countPhenotype[10]++;
            }
            else if(str.equals(myDataSet[i][4])  && myDataSet[i][7].equals("integral membrane")){
                countPhenotype[11]++;
            }
            else if(str.equals(myDataSet[i][4])  && myDataSet[i][7].equals("cell wall")){
                countPhenotype[12]++;
            }
            else if(str.equals(myDataSet[i][4])  && myDataSet[i][7].equals("lipid particles")){
                countPhenotype[13]++;
            }
            else if(str.equals(myDataSet[i][4])  && myDataSet[i][7].equals("transport vesicles")){
                countPhenotype[14]++;
            }
        }
    }

    private void findProbMotif(String str){

        setupFunction(countMotif);
        for(int i = 0; i < totalData; i++){
            if(str.equals(myDataSet[i][5])  && myDataSet[i][7].equals("nucleus")){
                countMotif[0]++;
            }
            else if(str.equals(myDataSet[i][5])  && myDataSet[i][7].equals("cytoplasm")){
                countMotif[1]++;
            }
            else if(str.equals(myDataSet[i][5])  && myDataSet[i][7].equals("cytoskeleton")){
                countMotif[2]++;
            }
            else if(str.equals(myDataSet[i][5])  && myDataSet[i][7].equals("mitochondria")){
                countMotif[3]++;
            }
            else if(str.equals(myDataSet[i][5])  && myDataSet[i][7].equals("plasma membrane")){
                countMotif[4]++;
            }
            else if(str.equals(myDataSet[i][5])  && myDataSet[i][7].equals("ER")){
                countMotif[5]++;
            }
            else if(str.equals(myDataSet[i][5])  && myDataSet[i][7].equals("golgi")){
                countMotif[6]++;
            }
            else if(str.equals(myDataSet[i][5])  && myDataSet[i][7].equals("vacuole")){
                countMotif[7]++;
            }
            else if(str.equals(myDataSet[i][5])  && myDataSet[i][7].equals("peroxisome")){
                countMotif[8]++;
            }
            else if(str.equals(myDataSet[i][5])  && myDataSet[i][7].equals("endosome")){
                countMotif[9]++;
            }
            else if(str.equals(myDataSet[i][5])  && myDataSet[i][7].equals("extracellular")){
                countMotif[10]++;
            }
            else if(str.equals(myDataSet[i][5])  && myDataSet[i][7].equals("integral membrane")){
                countMotif[11]++;
            }
            else if(str.equals(myDataSet[i][5])  && myDataSet[i][7].equals("cell wall")){
                countMotif[12]++;
            }
            else if(str.equals(myDataSet[i][5])  && myDataSet[i][7].equals("lipid particles")){
                countMotif[13]++;
            }
            else if(str.equals(myDataSet[i][5])  && myDataSet[i][7].equals("transport vesicles")){
                countMotif[14]++;
            }
        }
    }

    private void findProbChromosome(String str){

        setupFunction(countChromosome);
        for(int i = 0; i < totalData; i++){
            if(str.equals(myDataSet[i][6])  && myDataSet[i][7].equals("nucleus")){
                countChromosome[0]++;
            }
            else if(str.equals(myDataSet[i][6])  && myDataSet[i][7].equals("cytoplasm")){
                countChromosome[1]++;
            }
            else if(str.equals(myDataSet[i][6])  && myDataSet[i][7].equals("cytoskeleton")){
                countChromosome[2]++;
            }
            else if(str.equals(myDataSet[i][6])  && myDataSet[i][7].equals("mitochondria")){
                countChromosome[3]++;
            }
            else if(str.equals(myDataSet[i][6])  && myDataSet[i][7].equals("plasma membrane")){
                countChromosome[4]++;
            }
            else if(str.equals(myDataSet[i][6])  && myDataSet[i][7].equals("ER")){
                countChromosome[5]++;
            }
            else if(str.equals(myDataSet[i][6])  && myDataSet[i][7].equals("golgi")){
                countChromosome[6]++;
            }
            else if(str.equals(myDataSet[i][6])  && myDataSet[i][7].equals("vacuole")){
                countChromosome[7]++;
            }
            else if(str.equals(myDataSet[i][6])  && myDataSet[i][7].equals("peroxisome")){
                countChromosome[8]++;
            }
            else if(str.equals(myDataSet[i][6])  && myDataSet[i][7].equals("endosome")){
                countChromosome[9]++;
            }
            else if(str.equals(myDataSet[i][6])  && myDataSet[i][7].equals("extracellular")){
                countChromosome[10]++;
            }
            else if(str.equals(myDataSet[i][6])  && myDataSet[i][7].equals("integral membrane")){
                countChromosome[11]++;
            }
            else if(str.equals(myDataSet[i][6])  && myDataSet[i][7].equals("cell wall")){
                countChromosome[12]++;
            }
            else if(str.equals(myDataSet[i][6])  && myDataSet[i][7].equals("lipid particles")){
                countChromosome[13]++;
            }
            else if(str.equals(myDataSet[i][6])  && myDataSet[i][7].equals("transport vesicles")){
                countChromosome[14]++;
            }
        }
    }

    private void findProbGeneID(String[][] testData){


        for(int i = 0; i < testData.length;i++) {
            String name = "";
            if(!uniqueLocalization.contains(testData[i][0])) {

                for (int j = 0; j < testData.length; j++) {

                    if (testData[i][0].equals(testData[j][0])&&testData[j][1].equals("nucleus")) {
                        countGeneID[0]++;
                    } else if (testData[i][0].equals(testData[j][0])&&testData[j][1].equals("cytoplasm")) {
                        countGeneID[1]++;
                    } else if (testData[i][0].equals(testData[j][0])&&testData[j][1].equals("cytoskeleton")) {
                        countGeneID[2]++;
                    } else if (testData[i][0].equals(testData[j][0])&&testData[j][1].equals("mitochondria")) {
                        countGeneID[3]++;
                    } else if (testData[i][0].equals(testData[j][0])&&testData[j][1].equals("plasma membrane")) {
                        countGeneID[4]++;
                    } else if (testData[i][0].equals(testData[j][0])&&testData[j][1].equals("ER")) {
                        countGeneID[5]++;
                    } else if (testData[i][0].equals(testData[j][0])&&testData[j][1].equals("golgi")) {
                        countGeneID[6]++;
                    } else if (testData[i][0].equals(testData[j][0])&&testData[j][1].equals("vacuole")) {
                        countGeneID[7]++;
                    } else if (testData[i][0].equals(testData[j][0])&&testData[j][1].equals("peroxisome")) {
                        countGeneID[8]++;
                    } else if (testData[i][0].equals(testData[j][0])&&testData[j][1].equals("endosome")) {
                        countGeneID[9]++;
                    } else if (testData[i][0].equals(testData[j][0])&&testData[j][1].equals("extracellular")) {
                        countGeneID[10]++;
                    } else if (testData[i][0].equals(testData[j][0])&&testData[j][1].equals("integral membrane")) {
                        countGeneID[11]++;
                    } else if (testData[i][0].equals(testData[j][0])&&testData[j][1].equals("cell wall")) {
                        countGeneID[12]++;
                    } else if (testData[i][0].equals(testData[j][0])&&testData[j][1].equals("lipid particles")) {
                        countGeneID[13]++;
                    } else if (testData[i][0].equals(testData[j][0])&&testData[j][1].equals("transport vesicles")) {
                        countGeneID[14]++;
                    }
                }
                name = setLocalName();
                setupFunction(countGeneID);
            }
            uniqueLocalization.put(testData[i][0], name);
        }
    }

    private String setLocalName() {

        //setupFunction(probabilityResult);

//        probabilityResult[0] = (this.countGeneID[0]/nucleus) * p_nucleus;
//        probabilityResult[1] = (this.countGeneID[1]/cytoplasm) * p_cytoplasm;
//        probabilityResult[2] = (this.countGeneID[2]/cytoskeleton)* p_cytoskeleton;
//        probabilityResult[3] = (this.countGeneID[3]/mitochondria)* p_mitochondria;
//        probabilityResult[4] = (this.countGeneID[4]/plasma_membrane) * p_plasma_membrane;
//        probabilityResult[5] = (this.countGeneID[5]/ER) * p_ER;
//        probabilityResult[6] = (this.countGeneID[6]/golgi)*  p_golgi;
//        probabilityResult[7] = (this.countGeneID[7]/vacuole)* p_vacuole;
//        probabilityResult[8] = (this.countGeneID[8]/peroxisome) * p_peroxisome;
//        probabilityResult[9] = (this.countGeneID[9]/endosome) * p_endosome;
//        probabilityResult[10] = (this.countGeneID[10]/extracellular)* p_extracellular;
//        probabilityResult[11] = (this.countGeneID[11]/integral_membrane)*  p_integral_membrane;
//        probabilityResult[12] = (this.countGeneID[12]/cell_wall) * p_cell_wall;
//        probabilityResult[13] = (this.countGeneID[13]/lipid_particles) *p_lipid_particles;
//        probabilityResult[14] = (this.countGeneID[14]/transport_vesicles) * p_transport_vesicles;

        String name;
        double largest = 0;
        int num = 0;
        for(int k = 0; k < countGeneID.length;k++){
            if(countGeneID[k] > largest){
                largest = countGeneID[k];
                num = k + 1;
            }
        }

        switch (num){
            case 1:  name = "nucleus"; break;
            case 2: name =  "cytoplasm";break;
            case 3: name =  "cytoskeleton";break;
            case 4: name =  "mitochondria";break;
            case 5: name = "plasma membrane";break;
            case 6: name = "ER";break;
            case 7: name = "golgi";break;
            case 8: name = "vacuole";break;
            case 9: name = "peroxisome";break;
            case 10: name = "endosome";break;
            case 11: name = "extracellular";break;
            case 12: name = "integral membrane";break;
            case 13: name = "cell wall";break;
            case 14: name = "lipid particles";break;
            case 15: name = "transport vesicles";break;

            default: name = "cytoplasm";break;
        }
        return name;
    }


    public Hashtable<String,String> findTheFinalLocalization(String[][] testData){
            findProbGeneID(testData);
         return uniqueLocalization;
    }



    private int calculateFinalProb(double[] probabilityResult){
        int chosen = 0;
        setupFunction(probabilityResult);

        probabilityResult[0] = (this.countEssential[0]/nucleus) * (this.countClass[0]/nucleus) * (this.countComplex[0]/nucleus) * (this.countMotif[0]/nucleus)*(this.countPhenotype[0]/nucleus) * (this.countChromosome[0]/nucleus)* p_nucleus;
        probabilityResult[1] = (this.countEssential[1]/cytoplasm) * (this.countClass[1]/cytoplasm) * (this.countComplex[1]/cytoplasm) * (this.countMotif[1]/cytoplasm)*(this.countPhenotype[1]/cytoplasm) * (this.countChromosome[1]/cytoplasm) * p_cytoplasm;
        probabilityResult[2] = (this.countEssential[2]/cytoskeleton) * (this.countClass[2]/cytoskeleton) * (this.countComplex[2]/cytoskeleton) * (this.countMotif[2]/cytoskeleton)*(this.countPhenotype[2]/cytoskeleton) * (this.countChromosome[2]/cytoskeleton)* p_cytoskeleton;
        probabilityResult[3] = (this.countEssential[3]/mitochondria)* (this.countClass[2]/mitochondria) * (this.countComplex[2]/mitochondria) * (this.countMotif[2]/mitochondria)*(this.countPhenotype[2]/mitochondria) * (this.countChromosome[2]/mitochondria) * p_mitochondria;
        probabilityResult[4] = (this.countEssential[4]/plasma_membrane)* (this.countClass[2]/plasma_membrane) * (this.countComplex[2]/plasma_membrane) * (this.countMotif[2]/plasma_membrane)*(this.countPhenotype[2]/plasma_membrane) * (this.countChromosome[2]/plasma_membrane) * p_plasma_membrane;
        probabilityResult[5] = (this.countEssential[5]/ER)* (this.countClass[2]/ER) * (this.countComplex[2]/ER) * (this.countMotif[2]/ER)*(this.countPhenotype[2]/ER) * (this.countChromosome[2]/ER) * p_ER;
        probabilityResult[6] = (this.countEssential[6]/golgi)* (this.countClass[2]/golgi) * (this.countComplex[2]/golgi) * (this.countMotif[2]/golgi)*(this.countPhenotype[2]/golgi) * (this.countChromosome[2]/golgi) * p_golgi;
        probabilityResult[7] = (this.countEssential[7]/vacuole)* (this.countClass[2]/vacuole) * (this.countComplex[2]/vacuole) * (this.countMotif[2]/vacuole)*(this.countPhenotype[2]/vacuole) * (this.countChromosome[2]/vacuole) * p_vacuole;
        probabilityResult[8] = (this.countEssential[8]/peroxisome) * (this.countClass[2]/peroxisome) * (this.countComplex[2]/peroxisome) * (this.countMotif[2]/peroxisome)*(this.countPhenotype[2]/peroxisome) * (this.countChromosome[2]/peroxisome)* p_peroxisome;
        probabilityResult[9] = (this.countEssential[9]/endosome)* (this.countClass[2]/endosome) * (this.countComplex[2]/endosome) * (this.countMotif[2]/endosome)*(this.countPhenotype[2]/endosome) * (this.countChromosome[2]/endosome) * p_endosome;
        probabilityResult[10] = (this.countEssential[10]/extracellular)* (this.countClass[2]/extracellular) * (this.countComplex[2]/extracellular) * (this.countMotif[2]/extracellular)*(this.countPhenotype[2]/extracellular) * (this.countChromosome[2]/extracellular) * p_extracellular;
        probabilityResult[11] = (this.countEssential[11]/integral_membrane)* (this.countClass[2]/integral_membrane) * (this.countComplex[2]/integral_membrane) * (this.countMotif[2]/integral_membrane)*(this.countPhenotype[2]/integral_membrane) * (this.countChromosome[2]/integral_membrane) * p_integral_membrane;
        probabilityResult[12] = (this.countEssential[12]/cell_wall) *(this.countClass[2]/cell_wall) * (this.countComplex[2]/cell_wall) * (this.countMotif[2]/cell_wall)*(this.countPhenotype[2]/cell_wall) * (this.countChromosome[2]/cell_wall) * p_cell_wall;
        probabilityResult[13] = (this.countEssential[13]/lipid_particles) * (this.countClass[2]/lipid_particles) * (this.countComplex[2]/lipid_particles) * (this.countMotif[2]/lipid_particles)*(this.countPhenotype[2]/lipid_particles) * (this.countChromosome[2]/lipid_particles) * p_lipid_particles;
        probabilityResult[14] = (this.countEssential[14]/transport_vesicles)* (this.countClass[2]/transport_vesicles) * (this.countComplex[2]/transport_vesicles) * (this.countMotif[2]/transport_vesicles)*(this.countPhenotype[2]/transport_vesicles) * (this.countChromosome[2]/transport_vesicles) * p_transport_vesicles;

        double number2compare = 0;
        for (int i = 0; i < probabilityResult.length; i++){
            if(probabilityResult[i] > number2compare){
                number2compare = probabilityResult[i];
                chosen = i+1;
            }
        }
        return chosen;
    }

}
