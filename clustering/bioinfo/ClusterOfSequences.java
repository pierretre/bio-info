package bioinfo;

import java.util.ArrayList;
import java.util.StringJoiner;
import java.util.Locale;

public class ClusterOfSequences {
	
	private ArrayList<ClusterOfSequences> subClusters = new ArrayList<ClusterOfSequences>();
	private ArrayList<Sequence> elements = new ArrayList<Sequence>();

	public ClusterOfSequences(Sequence element) {
		this.elements.add(element);
	}
	
	public ClusterOfSequences(ArrayList<Sequence> eltList) {
		this.elements.addAll(eltList);
	}
	
	public ClusterOfSequences(ClusterOfSequences cluster1, ClusterOfSequences cluster2) {
		this.subClusters.add(cluster1);
		this.subClusters.add(cluster2);
		this.elements.addAll(cluster1.elements);
		this.elements.addAll(cluster2.elements);
	}
	
	public String getNewickIntermediate() {
		StringJoiner joiner = new StringJoiner(", ");

		if(this.subClusters.isEmpty()){
			if(this.elements.size() == 1) return this.elements.get(0).toString();
			this.elements.forEach(seq -> joiner.add(seq.toString()));
		}
		else
			this.subClusters.forEach(cluster -> joiner.add(cluster.getNewickIntermediate()));
		return "("+joiner.toString()+")";
	}
	
	public String getNewick() {
		return this.getNewickIntermediate() + ";";
	}
	
	private double depth() {
		if(!this.subClusters.isEmpty()) {
			double maxDepth = 0.;
			for(ClusterOfSequences cluster : this.subClusters) {
				double depth = cluster.depth();
				if(depth > maxDepth) maxDepth = depth;
			}
			return maxDepth + 1.;
		}
		return 0.;
	}

	public String getNewickIntermediateAligned(double alignRight) {
		StringJoiner joiner = new StringJoiner(", ");

		if(this.subClusters.isEmpty()){
			if(this.elements.size() == 1) 
				return this.elements.get(0).toString() + ((alignRight > 0)? ":" + String.valueOf(alignRight + 1) : "");
			this.elements.forEach(seq -> joiner.add(seq.toString() + ((alignRight > 0)? ":" + String.valueOf(alignRight + 1) : "")));
		}
		else
			this.subClusters.forEach(cluster -> joiner.add(cluster.getNewickIntermediateAligned(alignRight - 1)));
		return "("+joiner.toString()+")";
	}
	
	public String getNewickAligned() {
		return this.getNewickIntermediateAligned(this.depth()) + ";";
	}

	private double depthDistance(double distance) {
		if(!this.subClusters.isEmpty()) {
			double maxDepth = 0.;
			for(ClusterOfSequences cluster : this.subClusters) {
				double depth = cluster.depthDistance(1 - this.linkage(cluster));
				if(depth > maxDepth) maxDepth = depth;
			}
			return maxDepth + distance;
		}
		return 0.1;
	}

	public String getNewickIntermediateAlignedDistance(double alignRight, double distance) {
		StringJoiner joiner = new StringJoiner(", ");

		if(this.subClusters.isEmpty()){
			if(this.elements.size() == 1) 
				return this.elements.get(0).toString() + ":" + String.format(Locale.ENGLISH, "%.3f", alignRight + distance);
			this.elements.forEach(seq -> joiner.add(seq.toString() + ":" + String.format(Locale.ENGLISH, "%.3f", alignRight + distance)));
		}
		else
			this.subClusters.forEach(cluster -> {
				double cDist = 1 - this.linkage(cluster);
				joiner.add(cluster.getNewickIntermediateAlignedDistance(alignRight - cDist, cDist));
			});
		return "("+joiner.toString()+")" + ((distance != -1)? ":"+String.format(Locale.ENGLISH, "%.3f", distance) : "" );
	}
	
	public String getNewickAlignedDistance() {
		return this.getNewickIntermediateAlignedDistance(this.depthDistance(0), -1) + ";";
	}
	
	public double linkage(ClusterOfSequences aCluster) {
		double dist = 0.;
		
		for (Sequence seq1 : this.elements) {
			for (Sequence seq2 : aCluster.elements) {
				dist += 1 - seq1.distance(seq2);
			}
		}
		dist = dist / (this.elements.size() * aCluster.elements.size());
		return dist;
	}

	public void clusterize() {
		this.clusterizeAgglomerative();
		// this.clusterizeDivisive();
	}	

	public void clusterizeAgglomerative() {
		for (Sequence seq : this.elements) 
			this.subClusters.add(new ClusterOfSequences(seq));
		
		while(this.subClusters.size() > 2){
			double bestProximity = 0.;
			ClusterOfSequences closestCluster = null;
	
			for (ClusterOfSequences cluster1 : this.subClusters.subList(0, this.subClusters.size() - 1)) {
				for (ClusterOfSequences cluster2 : this.subClusters.subList(this.subClusters.indexOf(cluster1) + 1, this.subClusters.size())) {
					double linkage = cluster1.linkage(cluster2);
					if(linkage > bestProximity){
						bestProximity = linkage;
						closestCluster = new ClusterOfSequences(cluster1, cluster2);
					}					
				}
			}
			this.subClusters.removeAll(closestCluster.subClusters);
			this.subClusters.add(closestCluster);
		}
	}

	public void clusterizeDivisive() {		
		if(this.elements.size() == 2){
			this.subClusters.add(new ClusterOfSequences(this.elements.get(0)));
			this.subClusters.add(new ClusterOfSequences(this.elements.get(1)));
		}else if(this.elements.size() > 2){

			double maxDistance = 0.;
			Sequence farsestSeq1 = null;
			Sequence farsestSeq2 = null;
	
			for (Sequence seq1 : this.elements.subList(0, this.elements.size() - 1)) {
				for (Sequence seq2 : this.elements.subList(this.elements.indexOf(seq1) + 1, this.elements.size())) {
					double dist = seq1.distance(seq2);
					if(dist > maxDistance){
						maxDistance = dist;
						farsestSeq1 = seq1;
						farsestSeq2 = seq2;
					}					
				}
			}

			ArrayList<Sequence> otherElements = new ArrayList<>(this.elements);
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

			ClusterOfSequences cluster1 = new ClusterOfSequences(group1);
			ClusterOfSequences cluster2 = new ClusterOfSequences(group2);

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
//		listSeq.add(seq6);
//		listSeq.add(seq7);
		
		ClusterOfSequences bioCluster = new ClusterOfSequences(listSeq);
		System.out.println(bioCluster.getNewick());
		bioCluster.clusterizeAgglomerative();
		System.out.println(bioCluster.getNewick());

		ClusterOfSequences bioCluster2 = new ClusterOfSequences(listSeq);
		bioCluster2.clusterizeDivisive();
		System.out.println(bioCluster2.getNewick());
	}

	public String toString(){
		return this.elements.toString();
	}
}
