package bioinfo;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.LinkedList;

public class TimeComparison {

	public static void main(String[] args) {
		try {
			
			/*
			 * DNA
			 */
			LinkedList<String> samples = new LinkedList<String>();
			samples.add("Homo sapiens (HBA1)");
			samples.add("Pan paniscus (HBA1)");
			samples.add("Pan troglodytes (HBA1)");
			samples.add("Mus musculus (Hba_a1)");
			samples.add("Mus musculus (Hba_a2)");
			samples.add("Rattus norvegicus (Hba_a2-1)");
			samples.add("Rattus norvegicus (Hba_a2-2)");
			samples.add("Rattus norvegicus (Hba_a3)");
			samples.add("Felis catus (HBA1)");
			samples.add("Bos taurus (HBA)");
			samples.add("Danio rerio (hbaa1)");
			samples.add("Macaca mulatta (HBA2)");
			samples.add("Xenopus tropicalis (hba1)");
			int nbSamples = samples.size();
			
			ArrayList<Sequence> dataNucleotides = new ArrayList<Sequence>(nbSamples);
			String seqFilePath;
			for (String currentSample : samples) {
				seqFilePath = System.getProperty("user.dir") + "/data/" + Utils.getTaxon(currentSample).replace(' ', '_') + "_" + Utils.getGeneName(currentSample) + "_sequence.fa";
				dataNucleotides.add(new SequenceLabeled(new File(seqFilePath), Utils.getTaxon(currentSample) + " " + Utils.getGeneName(currentSample)));
			}
			
			
			// ClusterOfSequences :
			ClusterOfSequences clusterHemoglobinNucleotides = new ClusterOfSequences(dataNucleotides);
			long start = System.nanoTime();
			clusterHemoglobinNucleotides.clusterizeAgglomerative();
			
			System.out.println();
			System.out.println("ClusterOfSequence execution time (agglomerative) : "+String.valueOf((System.nanoTime() - start) * Math.pow(10, -9))+" seconds");
			System.out.println(clusterHemoglobinNucleotides.getNewick());
			System.out.println();

			clusterHemoglobinNucleotides = new ClusterOfSequences(dataNucleotides);
			start = System.nanoTime();
			clusterHemoglobinNucleotides.clusterizeDivisive();
			
			System.out.println();
			System.out.println("ClusterOfSequence execution time (divisive) : "+String.valueOf((System.nanoTime() - start) * Math.pow(10, -9))+" seconds");
			System.out.println(clusterHemoglobinNucleotides.getNewick());
			System.out.println();
			
			// ClusterOfSequencesBis :
			ClusterOfSequencesBis clusterHemoglobinNucleotidesBis = new ClusterOfSequencesBis(dataNucleotides);
			
			start = System.nanoTime();
			clusterHemoglobinNucleotidesBis.clusterizeAgglomerative();
			
			System.out.println();
			System.out.println("ClusterOfSequenceBis execution time (agglomerative) : "+String.valueOf((System.nanoTime() - start) * Math.pow(10, -9))+" seconds");
			System.out.println(clusterHemoglobinNucleotidesBis.getNewick());
			System.out.println();

			clusterHemoglobinNucleotidesBis = new ClusterOfSequencesBis(dataNucleotides);
			
			start = System.nanoTime();
			clusterHemoglobinNucleotidesBis.clusterizeDivisive();
			
			System.out.println();
			System.out.println("ClusterOfSequenceBis execution time (divisive) : "+String.valueOf((System.nanoTime() - start) * Math.pow(10, -9))+" seconds");
			System.out.println(clusterHemoglobinNucleotidesBis.getNewick());
			System.out.println();
						
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	

}
