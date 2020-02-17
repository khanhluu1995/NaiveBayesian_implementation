import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Scanner;
import java.util.Set;

public class Main {
    public static void main(String args[]) throws FileNotFoundException {
        File file = new File("Genes_relation.data");
        Scanner scanner = new Scanner(file);
        scanner.nextLine();
        double totalExample = 0;
        String[][] trainingData,testData;

        Hashtable<String, String> geneIDTable;
        Hashtable<String, Integer> localizationTable = new Hashtable<>();


        while (scanner.hasNext()) {
            totalExample++;
            String[] splitLine = scanner.nextLine().split(",");

//            addGeneIDToHashtable(geneIDTable, splitLine);
//            addEssentialToHashtable(essentialTable,splitLine);
//            addClassToHashtable(classTable,splitLine);
//            addComplexToHashtable(complexTable,splitLine);
//            addPhenotypeToHashtable(phenotypeTable,splitLine);
//            addChromosomeToHashtable(chromosomeTable,splitLine);
//            addMotifToHashtable(motifTable,splitLine);
            addLocalizationToHashtable(localizationTable,splitLine);

        }

        int amountOfAttr = 8;
        trainingData = new String[(int) totalExample][amountOfAttr];
        scanner = new Scanner(file);
        scanner.nextLine();
        int i = 0;
        while (scanner.hasNext()) {

            String[] splitLine = scanner.nextLine().split(",");
            for(int j = 0; j < amountOfAttr; j++){
                if(j == 7){
                    trainingData[i][j] = splitLine[splitLine.length-1];
                }
                else {
                    trainingData[i][j] = splitLine[j];
                }
            }
            i++;
        }

        //print2D(trainingData);
        
        int totalTestExample = 0;

        NaiveBayesian_Implementation naiveBayesian_implementation = new NaiveBayesian_Implementation(localizationTable,trainingData ,totalExample);

        file = new File("Genes_relation.test");
        scanner = new Scanner(file);
        scanner.nextLine();

        //Create an output file
        File writeFileName = new File("output.txt");
        try{
            if (writeFileName.createNewFile()) {
                System.out.println("File created: " + writeFileName.getName());
            } else {
                System.out.println("File already exists.");
            }
        } catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        };

        try {
            FileWriter fileWriter = new FileWriter(writeFileName);
            fileWriter.write("GeneID,Localization" + "\n");
            int line = 0;
            while(scanner.hasNext()) {
                scanner.nextLine();
                totalTestExample++;
            }

            testData = new String[totalTestExample][2];

            scanner = new Scanner(file);
            scanner.nextLine();
            while(scanner.hasNext()){
                String[] splitLine = scanner.nextLine().split(",");
                //System.out.println(line + ": " + naiveBayesian_implementation.predictLocalization(splitLine));
                testData[line][0]  = splitLine[0];
                testData[line][1]  = naiveBayesian_implementation.predictLocalization(splitLine);
                line++;
            }

            geneIDTable = naiveBayesian_implementation.findTheFinalLocalization(testData);
//            print2D(testData);
//            System.out.println(testData.length);
//            System.out.println(geneIDTable.toString());
//            System.out.println(geneIDTable.size());

            Set<String> keys = geneIDTable.keySet();
            String[] geneID = new String[keys.size()];
            String[] localization = new String[keys.size()];
            i = 0;
            for(String key : keys){
                geneID[i] = key;
                localization[i] = geneIDTable.get(key);
                i++;
            }

            for(i =0; i < geneID.length;i++){
                fileWriter.write(geneID[i] + "," + localization[i] + "\n");
            }

            compareOutputVsKeys(geneID, localization);

            fileWriter.close();
        }catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }
    }

    private static void compareOutputVsKeys(String[] geneID, String[] localization) throws FileNotFoundException {
        File file;
        Scanner scanner;
        int i;
        double correct = 0;
        file = new File("keys.txt");
        scanner = new Scanner(file);
        scanner.nextLine();
        while(scanner.hasNext()){
            String[] splitLine = scanner.nextLine().split(",");
            for(i = 0; i < geneID.length;i++){
                if(splitLine[0].equals(geneID[i])){
                    if(splitLine[1].equals(localization[i])){
                        correct++;
                    }
                }
            }
        }

        double accuracy = correct/geneID.length;
        System.out.println("accuracy is: " +Math.round(accuracy * 100)  + "%");
    }

    public static void print2D(String mat[][])
    {
        // Loop through all rows
        for (int i = 0; i < mat.length; i++) {
            System.out.println();
            // Loop through all elements of current row
            for (int j = 0; j < mat[i].length; j++) {
                System.out.print(mat[i][j] + " ");
            }
        }
    }



    private static void addLocalizationToHashtable(Hashtable<String, Integer> localizationTable, String[] splitLine) {
        int rightValue = splitLine.length-1;


        if (localizationTable.containsKey(splitLine[rightValue])) {
            localizationTable.put(splitLine[rightValue], localizationTable.get(splitLine[rightValue])+1);
        } else {
            localizationTable.put(splitLine[rightValue], 1);
        }
    }





}
