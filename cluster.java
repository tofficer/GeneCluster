//Bioinformatics
//Author: Tyler Officer

import java.io.*;
import java.util.*;
import java.math.*;

public class cluster {

	static class Gene implements Comparable<Gene> {
		String id;
		String description;
		ArrayList<Double> vals;
		double avg;

		public Gene(String id, String description, ArrayList<Double> vals, double avg) {
			this.id = id;
			this.description = description;
			this.vals = vals;
			this.avg = avg;
		}

		@Override
		public int compareTo(Gene o) {
			if (this.avg <= o.avg) return -1;
			else return 1;
		}
	}

	static class Cluster implements Comparable<Cluster> {
		ArrayList<Gene> genes; //set of genes in the cluster
		HashMap<Cluster, Double> map; //distances to other clusters
		double avg;

		public Cluster(ArrayList<Gene> genes, double avg) {
			this.genes = genes;
			this.avg = avg;
			this.map = new HashMap<>();
		}

		@Override
		public int compareTo(Cluster o) {
			if (this.avg <= o.avg) return -1;
			else return 1;
		}
	}

	public static void main(String[] args) throws IOException {
		String infile = args[0];
		String ch = args[1];
		int k = Integer.valueOf(args[2]);

		ArrayList<Cluster> clusters = parse(infile);
		gene_cluster(clusters, ch, k);
	}

	//agglomerative hierarchical clustering algorithm
	public static void gene_cluster(ArrayList<Cluster> clusters, String ch, int k) {		
		//setup original distance matrix based on euclidian distance
		clusters = get_distances(clusters);

		while (clusters.size() > k) {
			//select the best pair of clusters to merge, then merge them
			Cluster[] selected = select(clusters);			
			clusters = merge(clusters, selected, ch);
		}

		print(clusters);
	}

	public static ArrayList<Cluster> get_distances(ArrayList<Cluster> clusters) {
		int n = clusters.size();
		for (int i = 0; i < n; i++) {
			Cluster c1 = clusters.get(i);
			Gene g1 = c1.genes.get(0); //initially only one gene per cluster
			
			for (int j = i+1; j < n; j++) {
				Cluster c2 = clusters.get(j);
				Gene g2 = c2.genes.get(0);
				
				double x = 0;
				for (int m = 0; m < g1.vals.size(); m++) {
					x += Math.pow((g1.vals.get(m)-g2.vals.get(m)), 2);
				}

				c1.map.put(c2, Math.sqrt(x));
				c2.map.put(c1, Math.sqrt(x));
			}
		}

		return clusters;
	}

	public static Cluster[] select(ArrayList<Cluster> clusters) {
		double min = Double.MAX_VALUE;
		Cluster[] selected = new Cluster[2];
		
		int n = clusters.size();
		for (int i = 0; i < n; i++) {
			Cluster c1 = clusters.get(i);
			for (int j = i+1; j < n; j++) {
				Cluster c2 = clusters.get(j);
				if (c1.map.get(c2) < min) {
					min = c1.map.get(c2);
					selected[0] = c1;
					selected[1] = c2;
				}
			}
		}

		return selected;
	}

	public static ArrayList<Cluster> merge(ArrayList<Cluster> clusters, Cluster[] selected, String ch) {
		Cluster c1 = selected[0];
		Cluster c2 = selected[1];

		for (Cluster curr : clusters) {
			if ((curr == c1) || (curr == c2)) continue;
			
			double update = -1;
			switch (ch) {
				case "S": update = single_link(curr, c1, c2); break;
				case "C": update = complete_link(curr, c1, c2); break;
				case "A": update = average_link(curr, c1, c2); break;
			}

			curr.map.put(c1, update);			
			c1.map.put(curr, update);

			curr.map.remove(c2);
		}
		c1.map.remove(c2);

		c1.avg = ( (c1.avg * c1.genes.size()) + (c2.avg * c2.genes.size()) ) / (c1.genes.size() + c2.genes.size());
		c1.genes.addAll(c2.genes);

		clusters.remove(c2);

		return clusters;
	}

	public static double single_link(Cluster curr, Cluster c1, Cluster c2) {
		return Math.min(curr.map.get(c1), curr.map.get(c2));
	}

	public static double complete_link(Cluster curr, Cluster c1, Cluster c2) {
		return Math.max(curr.map.get(c1), curr.map.get(c2));
	}
	
	public static double average_link(Cluster curr, Cluster c1, Cluster c2) {
		return ( (c1.genes.size()*curr.map.get(c1)) + (c2.genes.size()*curr.map.get(c2)) ) / (c1.genes.size()+c2.genes.size());
	}

	public static void print(ArrayList<Cluster> clusters) {		
		Collections.sort(clusters);
		for (Cluster cluster : clusters) {
			
			Collections.sort(cluster.genes);
			for (Gene g : cluster.genes) {
				BigDecimal db1 = new BigDecimal(g.avg).setScale(3, BigDecimal.ROUND_HALF_EVEN);
				System.out.println(g.id + " " + g.description + " " + db1);
			}
			
			BigDecimal db2 = new BigDecimal(cluster.avg).setScale(3, BigDecimal.ROUND_HALF_EVEN);
			System.out.println(db2);
			System.out.println();
		}
	}

	public static ArrayList<Cluster> parse(String infile) throws IOException, FileNotFoundException {
		ArrayList<Cluster> clusters = new ArrayList<>(); 
		try (BufferedReader br = new BufferedReader(new FileReader(infile))) {
			String line;
			while ((line = br.readLine()) != null) {
				String[] data = line.split("\t", 0);
				
				String id = data[0].trim();
				String description = data[1].trim();
				
				ArrayList<Double> vals = new ArrayList<>();
				double avg = 0;
				for (int i = 2; i < data.length; i++) {
					double t = Double.valueOf(data[i]);
					vals.add(t);
					avg += t;
				}
				avg /= vals.size();

				ArrayList<Gene> genes = new ArrayList<>();
				genes.add(new Gene(id, description, vals, avg));
				clusters.add(new Cluster(genes, avg));
			}
		}

		return clusters;
	}
}