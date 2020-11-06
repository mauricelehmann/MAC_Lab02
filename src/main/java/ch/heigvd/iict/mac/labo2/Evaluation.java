package ch.heigvd.iict.mac.labo2;

import org.apache.lucene.analysis.Analyzer;
import org.apache.lucene.analysis.StopwordAnalyzerBase;
import org.apache.lucene.analysis.core.WhitespaceAnalyzer;
import org.apache.lucene.analysis.en.EnglishAnalyzer;
import org.apache.lucene.analysis.standard.StandardAnalyzer;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

public class Evaluation {

    private static Analyzer analyzer = null;

    private static void readFile(String filename, Function<String, Void> parseLine)
            throws IOException {
        try (BufferedReader br = new BufferedReader(
                new InputStreamReader(
                        new FileInputStream(filename),
                        StandardCharsets.UTF_8)
        )) {
            String line = br.readLine();
            while (line != null) {
                line = line.trim();
                if (!line.isEmpty()) {
                    parseLine.apply(line);
                }
                line = br.readLine();
            }
        }
    }

    /*
     * Reading CACM queries and creating a list of queries.
     */
    private static List<String> readingQueries() throws IOException {
        final String QUERY_SEPARATOR = "\t";

        List<String> queries = new ArrayList<>();

        readFile("evaluation/query.txt", line -> {
            String[] query = line.split(QUERY_SEPARATOR);
            queries.add(query[1]);
            return null;
        });
        return queries;
    }

    /*
     * Reading stopwords
     */
    private static List<String> readingCommonWords() throws IOException {
        List<String> commonWords = new ArrayList<>();

        readFile("common_words.txt", line -> {
            commonWords.add(line);
            return null;
        });
        return commonWords;
    }


    /*
     * Reading CACM qrels and creating a map that contains list of relevant
     * documents per query.
     */
    private static Map<Integer, List<Integer>> readingQrels() throws IOException {
        final String QREL_SEPARATOR = ";";
        final String DOC_SEPARATOR = ",";

        Map<Integer, List<Integer>> qrels = new HashMap<>();

        readFile("evaluation/qrels.txt", line -> {
            String[] qrel = line.split(QREL_SEPARATOR);
            int query = Integer.parseInt(qrel[0]);

            List<Integer> docs = qrels.get(query);
            if (docs == null) {
                docs = new ArrayList<>();
            }

            String[] docsArray = qrel[1].split(DOC_SEPARATOR);
            for (String doc : docsArray) {
                docs.add(Integer.parseInt(doc));
            }

            qrels.put(query, docs);
            return null;
        });
        return qrels;
    }

    public static void main(String[] args) throws IOException {

        ///
        /// Reading queries and queries relations files
        ///
        List<String> queries = readingQueries();
        System.out.println("Number of queries: " + queries.size());

        Map<Integer, List<Integer>> qrels = readingQrels();
        System.out.println("Number of qrels: " + qrels.size());

        double avgQrels = 0.0;
        for (int q : qrels.keySet()) {
            avgQrels += qrels.get(q).size();
        }
        avgQrels /= qrels.size();
        System.out.println("Average number of relevant docs per query: " + avgQrels);

        //TODO student: use this when doing the english analyzer + common words
        List<String> commonWords = readingCommonWords();

        ///
        ///  Part I - Select an analyzer
        ///
        // TODO student: compare Analyzers here i.e. change analyzer to
        // the asked analyzers once the metrics have been implemented
        analyzer = new WhitespaceAnalyzer();


        ///
        ///  Part I - Create the index
        ///
        Lab2Index lab2Index = new Lab2Index(analyzer);
        lab2Index.index("documents/cacm.txt");

        ///
        ///  Part II and III:
        ///  Execute the queries and assess the performance of the
        ///  selected analyzer using performance metrics like F-measure,
        ///  precision, recall,...
        ///

        // TODO student
        // compute the metrics asked in the instructions
        // you may want to call these methods to get:
        // -  The query results returned by Lucene i.e. computed/empirical
        //    documents retrieved
        //        List<Integer> retrivedDocuments = lab2Index.search(query);
        //
        // - The true query results from qrels file i.e. genuine documents
        //   returned matching a query
        //        List<Integer> qrelResults = qrels.get(queryNumber);

        int queryNumber = 0;
        int totalRelevantDocs = 0;
        int totalRetrievedDocs = 0;
        int totalRetrievedRelevantDocs = 0;
        double avgPrecision = 0.0;
        double avgRPrecision = 0.0;
        double avgRecall = 0.0;
        double meanAveragePrecision = 0.0;
        double fMeasure = 0.0;

        // average precision at the 11 recall levels (0,0.1,0.2,...,1) over all queries
        double[] avgPrecisionAtRecallLevels = createZeroedRecalls();

        int nbQueries = queries.size();
        //int nbQueries = 1;
        //For every queries :
        for (queryNumber = 0; (queryNumber < nbQueries); queryNumber++) {

            //Get the docs results from the search
            List<Integer> retrivedDocuments = lab2Index.search(queries.get(queryNumber));
            //Skip the empty query results
            if(retrivedDocuments.size() == 0){
                continue;
            }
            //Get the "true" revelants docs for the query
            List<Integer> qrel = qrels.get(queryNumber + 1);
            //Some queries do not have "true" revelants docs
            if(qrel == null) qrel = new ArrayList<>();

            //Sum the relevants docs
            totalRelevantDocs += qrel.size();
            //Sum the retrieved documents
            totalRetrievedDocs += retrivedDocuments.size();


            //Compute average precision for a given query

            // Will be, for this example of docs: [] [x] [] [] [x] (where [x] are revelant docs)
            // AP == [(1 / 2) + (2 / 5)] / 2
            // Recall == 1 if theire is 2 true revelants docs ((# of revelant documents founded in the retrieved document / # of "true" revelant documents)

            //We'll use an iterator for better time complexity
            Iterator<Integer> iterator = retrivedDocuments.iterator();
            int revelantDocsCounter = 0;
            double singleQueryAVGPrecision = 0.0;
            double singleQueryRecall = 0.0;
            double rPrecision = 0.0;
            ArrayList<Double> precisionArray = new ArrayList<>();
            ArrayList<Double> recallArray = new ArrayList<>();


            //Iterate on all retrieved documents
            //Stop when : we iterate all the retrieved documents OR when we found all the revelant docs
            for(int docIdx = 1; (docIdx <= retrivedDocuments.size()) || (revelantDocsCounter < qrel.size()); docIdx++){
                if(!iterator.hasNext()) break;
                //If the doc in the retDocs list is a revelant document
                if(qrel.contains(iterator.next())){
                    revelantDocsCounter++;
                    singleQueryAVGPrecision += (revelantDocsCounter / (double)docIdx);
                    //Calculate the R-precision numerator
                    if(docIdx <= qrel.size()){
                        rPrecision += 1;
                    }
                }
                double recall = revelantDocsCounter / (double) qrel.size();
                double precision = revelantDocsCounter / (double) docIdx;
                precisionArray.add(precision);
                recallArray.add(recall);
            }
            //We define that if their is no "true" revelants documents at all, recall & AP are zero
            if(qrel.size() == 0) {
                singleQueryAVGPrecision = 0;
                singleQueryRecall = 0;
                rPrecision = 0;
            } else {
                singleQueryAVGPrecision = singleQueryAVGPrecision / (double) qrel.size();
                singleQueryRecall = revelantDocsCounter / (double) qrel.size();
                rPrecision = rPrecision / (double) qrel.size();
            }

            //Compute the Precision At Recall Levels
            double[] precisionAtRecallLevels = getAvgPrecisionAtRecallLevels(precisionArray, recallArray);
            for (int level= 0; level <= 10; level++) {
                avgPrecisionAtRecallLevels[level] += precisionAtRecallLevels[level];
            }

            //Sum the AP in the MAP
            meanAveragePrecision += singleQueryAVGPrecision;

            //Sum the recalls
            avgRecall += singleQueryRecall;

            //Precision is number of revelant doc in the query / all the docs (eg. 6/10)
            avgPrecision += (revelantDocsCounter/ ((double) retrivedDocuments.size()));
            totalRetrievedRelevantDocs += revelantDocsCounter;

            //Sum the R-Precision
            avgRPrecision += rPrecision;
        }
        //End of loop on every queries...

        //Get the average precision over all queries
        avgPrecision = avgPrecision / (double) nbQueries;
        //Get the mean average precision (MAP)
        meanAveragePrecision = meanAveragePrecision / (double) nbQueries;
        //Get the average recall
        avgRecall = avgRecall / (double) nbQueries;

        //Get the F-measure of average precision and average recall
        fMeasure = (2 * avgPrecision * avgRecall ) / (avgPrecision + avgRecall);

        //Get the average R-Precision
        avgRPrecision = avgRPrecision / (double) nbQueries;

        //Compute the average Precision At Recall Levels
        for (int level= 0; level <= 10; level++) {
            avgPrecisionAtRecallLevels[level] = avgPrecisionAtRecallLevels[level] / (double) nbQueries;
        }

        ///
        ///  Part IV - Display the metrics
        ///

        displayMetrics(totalRetrievedDocs, totalRelevantDocs,
                totalRetrievedRelevantDocs, avgPrecision, avgRecall, fMeasure,
                meanAveragePrecision, avgRPrecision,
                avgPrecisionAtRecallLevels);
    }
    private static void displayMetrics(
            int totalRetrievedDocs,
            int totalRelevantDocs,
            int totalRetrievedRelevantDocs,
            double avgPrecision,
            double avgRecall,
            double fMeasure,
            double meanAveragePrecision,
            double avgRPrecision,
            double[] avgPrecisionAtRecallLevels
    ) {
        String analyzerName = analyzer.getClass().getSimpleName();
        if (analyzer instanceof StopwordAnalyzerBase) {
            analyzerName += " with set size " + ((StopwordAnalyzerBase) analyzer).getStopwordSet().size();
        }
        System.out.println(analyzerName);

        System.out.println("Number of retrieved documents: " + totalRetrievedDocs);
        System.out.println("Number of relevant documents: " + totalRelevantDocs);
        System.out.println("Number of relevant documents retrieved: " + totalRetrievedRelevantDocs);

        System.out.println("Average precision: " + avgPrecision);
        System.out.println("Average recall: " + avgRecall);

        System.out.println("F-measure: " + fMeasure);

        System.out.println("MAP: " + meanAveragePrecision);

        System.out.println("Average R-Precision: " + avgRPrecision);

        System.out.println("Average precision at recall levels: ");
        for (int i = 0; i < avgPrecisionAtRecallLevels.length; i++) {
            System.out.println(String.format("\t%s: %s", i, avgPrecisionAtRecallLevels[i]));
        }
    }

    private static double[] createZeroedRecalls() {
        double[] recalls = new double[11];
        Arrays.fill(recalls, 0.0);
        return recalls;
    }

    private static double[] getAvgPrecisionAtRecallLevels(ArrayList<Double> precisionArray, ArrayList<Double> recallArray) {

        double[] recalls = new double[11];

        int index = 0;
        for (int i = 0; i <= 10 ; i ++) {

            double level = i / 10.0;
            // Find the first element of recall array that match with our level
            // We do this to avoid to re-iterate over all the precisionArray for every levels
            while(index < recallArray.size() && recallArray.get(index) < level){
                index++;
            }

            //Find the max value in the precisionArray, in the bounds : [index, size of the array]
            double max = 0.0;
            for (int j = index; j < precisionArray.size(); j++) {
                if(precisionArray.get(j) > max)
                    max = precisionArray.get(j);
            }
            recalls[i] = max;
        }

        return recalls;
    }

    private static double getAveragePrecision(List<Integer> retDocs, List<Integer> revDocs) {
        /*
        //We'll use an iterator for better time complexity
        Iterator<Integer> iterator = retDocs.iterator();
        int revelantDocsCounter = 0;
        double precision = 0.0;

        //Iterate on all retrieved documents
        //Stop when : we iterate all the retrieved documents OR when we found all the revelant docs
        for(int docIdx = 1; (docIdx <= retDocs.size()) || (revelantDocsCounter < revDocs.size()); docIdx++){
            //If the doc in the retDocs list is a revelant document
            if(revDocs.contains(iterator.next())){
                revelantDocsCounter++;
                precision += (revelantDocsCounter / (double)docIdx);
            }
        }
        return precision / revelantDocsCounter;
        */
        return 0.0;
    }
}