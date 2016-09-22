package pt.uminho.sysbio.biomass;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;
//import org.biojava.nbio.core.sequence.features.AbstractFeature;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;




public class ReadFile {

	private static final String SEP = ";";
	private static final String GENE_PREFIX = ">";
	private static final Logger logger = LoggerFactory.getLogger(ReadFile.class);
	
	private static int getStartGeneIndex(String[] words) {
		int startGeneIndex = -1;
		for (int i = 1; i < words.length; i++) {
			if (words[i].startsWith(GENE_PREFIX)) {
				startGeneIndex = i+1;
				break;
			}
		}
		return startGeneIndex;
	}
	
	private static int getEndGeneIndex(String[] words) {
		int startEndIndex = -1;
		for (int i = 1; i < words.length; i++) {
			if (words[i].endsWith(" ")) {
				startEndIndex = i;
				break;
			}
		}
		return startEndIndex;
	}
	
	private static String getGeneLocus(String[] words) {
		int startTaxIndex = getStartGeneIndex(words);
		if (startTaxIndex < 0) {
			return null;
		}
		
		List<String> wordList = new ArrayList<> ();
		//
		for (int i = startTaxIndex; i < words.length; i++) {
			wordList.add(words[i]);
		}
		return StringUtils.join(wordList, ' ');
	}
	
	private static String getName(String[] words) {
		int startTaxIndex = getStartGeneIndex(words);
		if (startTaxIndex < 0) {
			return null;
		}
		
		
		List<String> wordList = new ArrayList<> ();
		//
		for (int i = 1; i < startTaxIndex; i++) {
			wordList.add(words[i]);
		}
		return StringUtils.join(wordList, ' ');
	}
	
	private static Map<String, Double> readCsv(File path) {
		Map<String, Double> geneData = new HashMap<> ();
		try {
			InputStream is = new FileInputStream(path);
			List<String> lines = IOUtils.readLines(is);
			
			for (int i = 1; i < lines.size(); i++) {
				String line = lines.get(i);
				String[] col = line.split(SEP);
				geneData.put(col[0], Double.parseDouble(col[1]));
			}
			IOUtils.closeQuietly(is);
		} catch (IOException e) {
			e.printStackTrace();
		}
		return geneData;
	}
	
	private static Map<String, Double> readCsvTest() {
		Map<String, Double> geneData = new HashMap<> ();
		geneData.put("b0001", 30.0);
		geneData.put("b0002", 50.0);
		geneData.put("b0003", 20.0);
		return geneData;
	}
	
	
	
	
	public static Map<String, Double> processaProteinas1 (String AAseq, String GNdata)throws Exception {
		File file= new File (AAseq); // ir buscar caminho para o ficheiro � caixa DNAseq
		File filegene= new File (GNdata); // ir buscar caminho para o ficheiro � caixa DNAseq
		int tamanho=0;//n total de amino�cidos
		Object [] ret=new Object[2];

		Map<String, ProteinSequence> sequencias=FastaReaderHelper.readFastaProteinSequence(file); //ler ficheiro Fasta
		Map<String, Double> geneData = readCsv(new File(GNdata));
//		Map<String, Double> geneData = readCsvTest();
		Map<String, Map<String, Integer>> h = new HashMap<>(); //Hash com o AA e o numero de cada amino�cidos
		Hashtable<String, ProteinSequence> gene_seq = new Hashtable<> ();
		Map<String, ?> proteinsMap = new HashMap<> ();
		
		Map<String, Double> total = new HashMap<> ();
		
		for (String header : sequencias.keySet()) {
			String[] words = header.split("\\s+");
			String name = getName(words);
			String geneLocus = words [0];
			
			ProteinSequence proteinSequence = sequencias.get(header);

			//proteinsMap � uma tabela que contem o locus_tag e a respectiva sequencia de amino�cidos 

			if (geneLocus != null) {
				proteinsMap.put(geneLocus, null);
				System.out.print ("\n" + geneLocus);
				System.out.print("\t" + proteinSequence);
				
				gene_seq.put(geneLocus,proteinSequence);
			}
			else {
				logger.error("No locus tag for {}", name);
			}
			
			List<AminoAcidCompound> AAcomp = proteinSequence.getCompoundSet().getAllCompounds();
			Map<String, Integer> ocurrencias = new HashMap<> ();
			Map<String, Double> expressionMap = new HashMap<> ();
			
			for (AminoAcidCompound aa : AAcomp) {
				double expression = geneData.get(geneLocus);
				String base = aa.getBase();
				int aacount = StringUtils.countMatches(proteinSequence.getSequenceAsString(), aa.getBase());
				ocurrencias.put(base, aacount);
				expressionMap.put(base, aacount * expression);
//				System.out.println(aa.getBase() + " -> " + ocurrencias.get(aa.getBase()));
				if (!total.containsKey(base)) {
					total.put(base, 0.0);
				}
				double aaTotal = total.get(base) + (aacount * expression);
				total.put(base, aaTotal);
			}
			
			h.put(geneLocus, ocurrencias);
			
			int tamanhoseq=proteinSequence.getLength();
			System.out.print("\n" + tamanhoseq + "\n");
		
		
		
		
		}
		
		double t = 0;
		for (Double v : total.values()) {
			t+=v;
		}
		
		Map<String, Double> net = new HashMap<> ();
		for (String base : total.keySet()) {
			double v = total.get(base);
			net.put(base, v / t);
		}
		
		return net;
	}

	

	
	public static void main(String[] args) {
		try {
			processaProteinas1("C:/Users/Sophia/Desktop/Ferramenta/Cyanothece.faa","C:/Users/Sophia/Desktop/Ferramenta/data_cya_light.csv");
//			processaProteinas1("C:/Users/Sophia/Desktop/Ferramenta/seq_prot_test.faa","C:/Users/Sophia/Desktop/Ferramenta/data.csv");
			Map<String, String> geneT_faa = new HashMap<> ();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static void processaProteinas(String string, String string2) {
		// TODO Auto-generated method stub
		
	}
}

