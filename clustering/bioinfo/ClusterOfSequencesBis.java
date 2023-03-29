package bioinfo;

import java.util.ArrayList;
import java.util.StringJoiner;

public class ClusterOfSequencesBis {
	
	private ArrayList<ClusterOfSequencesBis> subClusters = new ArrayList<ClusterOfSequencesBis>();
	
	// added by me :
	private Sequence sequence = null;
	
	public ClusterOfSequencesBis(Sequence element) {
		this.sequence = element;
	}
	
	public ClusterOfSequencesBis(ArrayList<Sequence> eltList) {		
		for (Sequence seq : eltList) {
			this.subClusters.add(new ClusterOfSequencesBis(seq));
		}
	}
	
	public ClusterOfSequencesBis(ClusterOfSequencesBis cluster1, ClusterOfSequencesBis cluster2) {
		this.subClusters.add(cluster1);
		this.subClusters.add(cluster2);
	}
	
	public String getNewickIntermediate() {
		String newick = "";

		if(this.sequence != null)
			newick = this.sequence.toString();
		else{
			StringJoiner joiner = new StringJoiner(", ");
			this.subClusters.forEach(cluster -> joiner.add(cluster.getNewickIntermediate()));
			newick = "("+joiner.toString()+")";
		}
		return newick;
	}
	
	public String getNewick() {
		return this.getNewickIntermediate() + ";";
	}
	
	private ArrayList<Sequence> getAllElements(){
		ArrayList<Sequence> ret = new ArrayList<Sequence>();
		if(this.subClusters.isEmpty()) 
			ret.add(this.sequence);
		else {
			for(ClusterOfSequencesBis cluster : this.subClusters) 
				ret.addAll(cluster.getAllElements());
		}
		return ret;
	}
	
	public double linkage(ClusterOfSequencesBis aCluster) {
		double dist = 0.;
		
		ArrayList<Sequence> cElements = this.getAllElements();
		ArrayList<Sequence> aElements = aCluster.getAllElements();
		
		for (Sequence seq1 : cElements) {
			for (Sequence seq2 : aElements) {
				dist += 1 - seq1.distance(seq2);
			}
		}
		dist = dist / (cElements.size() * aElements.size());
		return dist;
	}
	
	public void clusterize() {
		this.clusterizeAgglomerative();
		// this.clusterizeDivisive();
	}	

	public void clusterizeAgglomerative() {
		
		while(this.subClusters.size() > 2){
			double bestProximity = 0.;
			ClusterOfSequencesBis closestCluster = null;
	
			for (ClusterOfSequencesBis cluster1 : this.subClusters.subList(0, this.subClusters.size() - 1)) {
				for (ClusterOfSequencesBis cluster2 : this.subClusters.subList(this.subClusters.indexOf(cluster1) + 1, this.subClusters.size())) {
					double linkage = cluster1.linkage(cluster2);
					if(linkage > bestProximity){
						bestProximity = linkage;
						closestCluster = new ClusterOfSequencesBis(cluster1, cluster2);
					}					
				}
			}
			this.subClusters.removeAll(closestCluster.subClusters);
			this.subClusters.add(closestCluster);
		}
	}

	public void clusterizeDivisive() {	
		ArrayList<Sequence> cElements = this.getAllElements();

		if(cElements.size() == 2){
			this.subClusters.add(new ClusterOfSequencesBis(cElements.get(0)));
			this.subClusters.add(new ClusterOfSequencesBis(cElements.get(1)));
		}else if(cElements.size() > 2){

			double maxDistance = 0.;
			Sequence farsestSeq1 = null;
			Sequence farsestSeq2 = null;
	
			for (Sequence seq1 : cElements.subList(0, cElements.size() - 1)) {
				for (Sequence seq2 : cElements.subList(cElements.indexOf(seq1) + 1, cElements.size())) {
					double dist = seq1.distance(seq2);
					if(dist > maxDistance){
						maxDistance = dist;
						farsestSeq1 = seq1;
						farsestSeq2 = seq2;
					}					
				}
			}

			ArrayList<Sequence> otherElements = new ArrayList<>(cElements);
			otherElements.remove(farsestSeq1);
			otherElements.remove(farsestSeq2);

			ArrayList<Sequence> group1 = new ArrayList<>();
			group1.add(farsestSeq1);

			ArrayList<Sequence> group2 = new ArrayList<>();
			group2.add(farsestSeq2);

			for (Sequence seq : otherElements) {
				if(seq.distance(farsestSeq1) < seq.distance(farsestSeq2))
					group1.add(seq);
				else
					group2.add(seq);
			}

			ClusterOfSequencesBis cluster1 = new ClusterOfSequencesBis(group1);
			ClusterOfSequencesBis cluster2 = new ClusterOfSequencesBis(group2);

			cluster1.clusterizeDivisive();
			cluster2.clusterizeDivisive();

			this.subClusters.add(cluster1);
			this.subClusters.add(cluster2);

		}
	}

	public static void main(String[] args) {
		Sequence seq1 = new Sequence("ATTACG");
		Sequence seq2 = new Sequence("ATATCG");
		Sequence seq3 = new Sequence("ACCCCG");
		Sequence seq4 = new Sequence("GCCGAG");
		Sequence seq5 = new Sequence("TCCCCG");
		Sequence seq6 = new Sequence("ATTAC");
		Sequence seq7 = new Sequence("ATATC");
		
		ArrayList<Sequence> listSeq = new ArrayList<Sequence>();
		listSeq.add(seq1);
		listSeq.add(seq2);
		listSeq.add(seq3);
		listSeq.add(seq4);
		listSeq.add(seq5);
		// listSeq.add(seq6);
		// listSeq.add(seq7);
		
		ClusterOfSequencesBis bioCluster = new ClusterOfSequencesBis(listSeq);
		System.out.println(bioCluster.getNewick());
		bioCluster.clusterize();
		System.out.println(bioCluster.getNewick());
	}
}
