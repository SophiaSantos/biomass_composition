package pt.uminho.sysbio.biomass;

import java.awt.EventQueue;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Spliterator;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;

import org.apache.commons.io.FileUtils;
import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.FastaReader;
import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.biojava3.core.sequence.io.FastaSequenceParser;
import org.eclipse.swt.SWT;
import org.eclipse.swt.widgets.Display;
import org.eclipse.swt.widgets.FileDialog;
import org.eclipse.swt.widgets.Shell;
import org.eclipse.swt.widgets.Text;


public class NucAaComp {

	private JFrame frame;
	private JLabel lblA_1;

	Display d= new Display();
	Shell s = new Shell(d);

	String DNAseq;
	String AAseq;
	String mRNAseq;
	String rRNAseq;
	String tRNAseq;
	double perc_mRNA;
	double perc_rRNA;
	double perc_tRNA;

	double bm_dna;
	double bm_rna;
	private JTabbedPane tabbedPane;
	private JPanel panel;
	private JPanel panel_1;
	private JPanel panel_2;
	private JTextField inputFP;
	private JTextField inputGD;
	private JTextField BM_prot;

	private JTextField inputFD;
	private JTextField BM_dna;
	private JTextField A_DNA;
	private JTextField T_DNA;
	private JTextField C_DNA;
	private JTextField G_DNA;
	private JTextField inputFmRNA;
	private JTextField inputFtRNA;
	private JTextField inputFrRNA;
	private JTextField BM_rna;
	private JTextField RNA_A;
	private JTextField RNA_U;
	private JTextField RNA_C;
	private JTextField RNA_G;
	Map<String, JTextField> aaMolTextFieldMap = new HashMap<> ();

	private JTextField dATP_mol;
	private JTextField dTTP_mol;
	private JTextField dCTP;
	private JTextField dGTP;
	private JTextField ATP_mol;
	private JTextField UTP_mol;
	private JTextField CTP_mol;
	private JTextField GTP_mol;
	private JTextField mRNA_per;
	private JTextField tRNA_per;
	private JTextField rRNA_per;

	final Text fileName = new Text(s, SWT.BORDER);
	
	private String[] aaArray = {
			"A", "Ala",
			"R", "Arg",
			"N", "Asn",
			"D", "Asp",
			"C", "Cys",
			"Q", "Gln",
			"E", "Glu",
			"G", "Gly",
			"H", "His",
			"I", "Ile",
			"L", "Leu",
			"K", "Lys",
			"M", "Met",
			"F", "Phe",
			"P", "Pro",
			"S", "Ser",
			"T", "Thr",
			"W", "Trp",
			"Y", "Tyr",
			"V", "Val",
	};
	
	private String[] dnaArray = {
			"A", "dATP",
			"C", "dCTP",
			"G", "dGTP",
			"T", "dTTP",
	};
	
	private String[] rnaArray = {
			"A", "ATP",
			"C", "CTP",
			"G", "GTP",
			"T", "UTP",
	};
	private HashMap<Object, Object> aagDWTextFieldMap;
	private Map<String, JTextField> aagDWMolTextFieldMap;
	private Map<String, JTextField> aaTextFieldMap;
	

	
	Map<String, HashMap <Double, Double>> protein = new HashMap<>();
	HashMap <Double, Double> pValues = new HashMap <>();
	
	Hashtable <String,Double> pContentm_m=new Hashtable <>();
	Hashtable <String,Double> pContentg_m=new Hashtable <>();
	Hashtable <String,Double> pContentm=new Hashtable <>();
	Hashtable<String, Double> pContentmmol_gDW = new Hashtable<String, Double> ();
	
	
	Hashtable <String,Double> dnaContentm_m=new Hashtable <>();
	Hashtable <String,Double> dnaContentg_m=new Hashtable <>();
	Hashtable <String,Double> dnaContentm=new Hashtable <>();
	Hashtable<String, Double> dnaContentmmol_gDW = new Hashtable<String, Double> ();

	
	Hashtable <String,Double> rnaContentm_m = new Hashtable <String, Double>();
	Hashtable <String,Double> rnaContentm_mr=new Hashtable <String, Double>();
	Hashtable <String,Double> rnaContentm_mm=new Hashtable <String, Double>();
	Hashtable <String,Double> rnaContentm_mt=new Hashtable <String, Double>();
	Hashtable<String, Double> rnaContentmmol_gDW = new Hashtable<String, Double> ();
	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					NucAaComp window = new NucAaComp();
					window.frame.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the application.
	 */
	public NucAaComp() {
		initialize();
	}

	/**
	 * Initialize the contents of the frame.
	 */
	private void initialize() {
		frame = new JFrame();
		frame.setBounds(100, 100, 800, 509);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.getContentPane().setLayout(null);

		tabbedPane = new JTabbedPane(JTabbedPane.TOP);
		tabbedPane.setBounds(0, 0, 774, 570);
		frame.getContentPane().add(tabbedPane);

		panel = new JPanel();
		tabbedPane.addTab("Protein", null, panel, null);
		panel.setLayout(null);

		JLabel lblProteinFastaFile = new JLabel("Protein FASTA file");
		lblProteinFastaFile.setBounds(10, 31, 123, 14);
		panel.add(lblProteinFastaFile);

		inputFP = new JTextField();
		inputFP.setColumns(10);
		inputFP.setBounds(128, 28, 113, 20);
		panel.add(inputFP);
		
		inputGD = new JTextField();
		inputGD.setColumns(10);
		inputGD.setBounds(128, 68, 113, 20);
		panel.add(inputGD);

		JButton OpenFP = new JButton("Open");
		OpenFP.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {

				FileDialog fd = new FileDialog(s, SWT.OPEN);
				fd.setText("Abrir");
				fd.setFilterPath("C:/");
				String[] filterExt = { "*.txt", "*.doc", ".rtf", "*.*" };
				fd.setFilterExtensions(filterExt);
				String selected = fd.open();
				inputFP.setText(selected);

			}
		});
		OpenFP.setBounds(251, 27, 72, 23);
		panel.add(OpenFP);

		BM_prot = new JTextField();
		BM_prot.setColumns(10);
		BM_prot.setBounds(673, 28, 86, 20);
		panel.add(BM_prot);

		JLabel lblProteinComposition = new JLabel("Cellular content Protein");
		lblProteinComposition.setBounds(518, 31, 145, 14);
		panel.add(lblProteinComposition);

		JButton DetermineP = new JButton("Determine");
		DetermineP.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {

				AAseq = inputFP.getText();
				//TODO: need work !!
				String geneData = inputGD.getText();
				double bm_prot = Double.parseDouble(BM_prot.getText());
				Map<String, Double> aaGeneData = null;
				try {
					System.out.println(AAseq + " " + geneData);
					Object[] res = processaProteinas();
					if (!geneData.trim().isEmpty()) {
						aaGeneData = ReadFile.processaProteinas1(AAseq, geneData);
					} else {
						
					}
					int tamanho=(Integer)res[1];
					@SuppressWarnings("unchecked")
					Hashtable <AminoAcidCompound, Integer> h1=(Hashtable<AminoAcidCompound, Integer>) res[0];

					Hashtable <String, Float> MWAA=new Hashtable <String, Float>();
					MWAA.put("A",  (float) 89.09);
					MWAA.put("R",  (float) 175.21);
					MWAA.put("N",  (float) 132.12);
					MWAA.put("D",  (float) 132.10);
					MWAA.put("C",  (float) 121.16);
					MWAA.put("Q",  (float) 146.15);
					MWAA.put("E",  (float) 146.12);
					MWAA.put("G",  (float) 75.07);
					MWAA.put("H",  (float) 155.16);
					MWAA.put("I",  (float) 131.18);
					MWAA.put("L",  (float) 131.18);
					MWAA.put("K",  (float) 147.20);
					MWAA.put("M",  (float) 149.21);
					MWAA.put("F",  (float) 165.19);
					MWAA.put("P",  (float) 115.13);
					MWAA.put("S",  (float) 105.09);
					MWAA.put("T",  (float) 119.12);
					MWAA.put("W",  (float) 204.23);
					MWAA.put("Y",  (float) 181.19);
					MWAA.put("V",  (float) 117.15); 

					Enumeration<AminoAcidCompound> em=h1.keys();
	//				System.out.println(h1);
					double soma=0;
					while (em.hasMoreElements()) {
						
						AminoAcidCompound key1=em.nextElement();

						if(h1.get(key1)!=0 && key1.getShortName()!="U" && key1.getShortName()!="X" && key1.getShortName()!="*"){
							
							System.out.println(key1.getShortName()+"-->"+pContentg_m.get(key1.getShortName()));
							
							if (!geneData.trim().isEmpty()) {
								pContentm.put(key1.getShortName(), aaGeneData.get(key1.getShortName()));
								double contentg=aaGeneData.get(key1.getShortName())*MWAA.get(key1.getShortName());
								pContentg_m.put(key1.getShortName(), contentg);
								soma+=contentg;
							}
							
							else {
								
							int count=h1.get(key1);
							double contentm_m=(float)count/(float)tamanho;
							pContentm.put(key1.getShortName(), contentm_m);
							double contentg=(float)contentm_m*MWAA.get(key1.getShortName());
							pContentg_m.put(key1.getShortName(), contentg);
							soma+=contentg;
							
							}

						}
					}

					double contentg_gA=pContentg_m.get("A")/soma;
					double contentg_gR=pContentg_m.get("R")/soma;
					double contentg_gN=pContentg_m.get("N")/soma;
					double contentg_gD=pContentg_m.get("D")/soma;
					double contentg_gC=pContentg_m.get("C")/soma;
					double contentg_gQ=pContentg_m.get("Q")/soma;
					double contentg_gE=pContentg_m.get("E")/soma;
					double contentg_gG=pContentg_m.get("G")/soma;
					double contentg_gH=pContentg_m.get("H")/soma;
					double contentg_gI=pContentg_m.get("I")/soma;
					double contentg_gL=pContentg_m.get("L")/soma;
					double contentg_gM=pContentg_m.get("M")/soma;
					double contentg_gK=pContentg_m.get("K")/soma;
					double contentg_gF=pContentg_m.get("F")/soma;
					double contentg_gP=pContentg_m.get("P")/soma;
					double contentg_gS=pContentg_m.get("S")/soma;
					double contentg_gT=pContentg_m.get("T")/soma;
					double contentg_gW=pContentg_m.get("W")/soma;
					double contentg_gY=pContentg_m.get("Y")/soma;
					double contentg_gV=pContentg_m.get("V")/soma;
					
					
					double contentmmol_gDW_A=(float)(contentg_gA*bm_prot*1000)/MWAA.get("A").shortValue();
					double contentmmol_gDW_R=(float)(contentg_gR*bm_prot*1000)/MWAA.get("R").shortValue();
					double contentmmol_gDW_N=(float)(contentg_gN*bm_prot*1000)/MWAA.get("N").shortValue();
					double contentmmol_gDW_D=(float)(contentg_gD*bm_prot*1000)/MWAA.get("D").shortValue();
					double contentmmol_gDW_C=(float)(contentg_gC*bm_prot*1000)/MWAA.get("C").shortValue();
					double contentmmol_gDW_Q=(float)(contentg_gQ*bm_prot*1000)/MWAA.get("Q").shortValue();
					double contentmmol_gDW_E=(float)(contentg_gE*bm_prot*1000)/MWAA.get("E").shortValue();
					double contentmmol_gDW_G=(float)(contentg_gG*bm_prot*1000)/MWAA.get("G").shortValue();
					double contentmmol_gDW_H=(float)(contentg_gH*bm_prot*1000)/MWAA.get("H").shortValue();
					double contentmmol_gDW_I=(float)(contentg_gI*bm_prot*1000)/MWAA.get("I").shortValue();
					double contentmmol_gDW_L=(float)(contentg_gL*bm_prot*1000)/MWAA.get("L").shortValue();
					double contentmmol_gDW_M=(float)(contentg_gM*bm_prot*1000)/MWAA.get("M").shortValue();
					double contentmmol_gDW_K=(float)(contentg_gK*bm_prot*1000)/MWAA.get("K").shortValue();
					double contentmmol_gDW_F=(float)(contentg_gF*bm_prot*1000)/MWAA.get("F").shortValue();
					double contentmmol_gDW_P=(float)(contentg_gP*bm_prot*1000)/MWAA.get("P").shortValue();
					double contentmmol_gDW_S=(float)(contentg_gS*bm_prot*1000)/MWAA.get("S").shortValue();
					double contentmmol_gDW_T=(float)(contentg_gT*bm_prot*1000)/MWAA.get("T").shortValue();
					double contentmmol_gDW_W=(float)(contentg_gW*bm_prot*1000)/MWAA.get("W").shortValue();
					double contentmmol_gDW_Y=(float)(contentg_gY*bm_prot*1000)/MWAA.get("Y").shortValue();
					double contentmmol_gDW_V=(float)(contentg_gV*bm_prot*1000)/MWAA.get("V").shortValue();

					pContentmmol_gDW.put("A", contentmmol_gDW_A);
					pContentmmol_gDW.put("R", contentmmol_gDW_R);
					pContentmmol_gDW.put("N", contentmmol_gDW_N);
					pContentmmol_gDW.put("D", contentmmol_gDW_D);
					pContentmmol_gDW.put("C", contentmmol_gDW_C);
					pContentmmol_gDW.put("Q", contentmmol_gDW_Q);
					pContentmmol_gDW.put("E", contentmmol_gDW_E);
					pContentmmol_gDW.put("G", contentmmol_gDW_G);
					pContentmmol_gDW.put("H", contentmmol_gDW_H);
					pContentmmol_gDW.put("I", contentmmol_gDW_I);
					pContentmmol_gDW.put("L", contentmmol_gDW_L);
					pContentmmol_gDW.put("M", contentmmol_gDW_M);
					pContentmmol_gDW.put("K", contentmmol_gDW_K);
					pContentmmol_gDW.put("F", contentmmol_gDW_F);
					pContentmmol_gDW.put("P", contentmmol_gDW_P);
					pContentmmol_gDW.put("S", contentmmol_gDW_S);
					pContentmmol_gDW.put("T", contentmmol_gDW_T);
					pContentmmol_gDW.put("W", contentmmol_gDW_W);
					pContentmmol_gDW.put("Y", contentmmol_gDW_Y);
					pContentmmol_gDW.put("V", contentmmol_gDW_V);
					
					for (String base : aaMolTextFieldMap.keySet()) {
						JTextField textField = aaMolTextFieldMap.get(base);
						textField.setText(String.format("%.5f", pContentm.get(base)));
					}
					
					for (String base : aaTextFieldMap.keySet()) {
						JTextField textField = aaTextFieldMap.get(base);
						textField.setText(String.format("%.5f", pContentmmol_gDW.get(base))); 
					}
					
					
					if (!geneData.trim().isEmpty()) {
						for (String base : aaMolTextFieldMap.keySet()) {
							JTextField textField = aaMolTextFieldMap.get(base);
							textField.setText(String.format("%.5f", aaGeneData.get(base)));
						}
					}
					
				} catch (Exception e1) {
					e1.printStackTrace();
				}
			}
		});
		DetermineP.setBounds(658, 59, 101, 23);
		panel.add(DetermineP);
		
		JPanel panel_3 = new JPanel();
		panel_3.setBounds(423, 162, 153, 262);
		panel.add(panel_3);
		
		JPanel panel_4 = new JPanel();
		panel_4.setBounds(586, 162, 150, 262);
		panel.add(panel_4);

		
		// colocar as caixas de texto dos AA (mmol/gDW)
		
				aaTextFieldMap = new HashMap<> ();
				for (int i = 0; i < aaArray.length / 2; i+=2) {
					String base = aaArray[i];
					String name = aaArray[i + 1];
					JTextField aaTextField = new JTextField();
					aaTextField.setColumns(10);
					
					JLabel aaLabel = new JLabel(name);
					
					panel_3.add(aaTextField);
					panel_3.add(aaLabel);
					
					aaTextFieldMap.put(base, aaTextField);
				}
				for (int i = aaArray.length / 2; i < aaArray.length; i+=2) {
					String base = aaArray[i];
					String name = aaArray[i + 1];
					JTextField aaTextField = new JTextField();
					aaTextField.setColumns(10);
					
					JLabel aaLabel = new JLabel(name);
					
					panel_4.add(aaTextField);
					panel_4.add(aaLabel);
					
					aaTextFieldMap.put(base, aaTextField);
				}


		JLabel lblBiomassCompositonAmino = new JLabel("Biomass Compositon Amino Acids (mmol/gDW)");
		lblBiomassCompositonAmino.setBounds(465, 132, 271, 14);
		panel.add(lblBiomassCompositonAmino);

		JLabel lblAminoAcidsContent = new JLabel("Amino Acids Content (mol/mol)");
		lblAminoAcidsContent.setBounds(112, 132, 180, 14);
		panel.add(lblAminoAcidsContent);

		JPanel panel_7 = new JPanel();
		panel_7.setBounds(198, 162, 150, 262);
		panel.add(panel_7);
		JPanel panel_8 = new JPanel();
		panel_8.setBounds(35, 162, 153, 262);
		panel.add(panel_8);
		
		JButton OpenGD = new JButton("Open");
		OpenGD.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				
				FileDialog gd = new FileDialog(s, SWT.OPEN);
				gd.setText("Abrir");
				gd.setFilterPath("C:/");
				String[] filterExt = { "*.txt", "*.doc", ".rtf", "*.*" };
				gd.setFilterExtensions(filterExt);
				String selected = gd.open();
				inputGD.setText(selected);
			}
		});
		
		JButton ExportP = new JButton("Export");
		ExportP.setText("Export");
		ExportP.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				FileDialog gd = new FileDialog(s, SWT.SAVE);
				gd.setText("Guardar");
				gd.setFilterPath("C:/");
				String[] filterExt = { "*.csv"};
				gd.setFilterExtensions(filterExt);
		        String selected = gd.open();
		        gd.setFileName(selected);
		        fileName.setText(selected);
		        FileWriter writer;
		        try {
					writer = new FileWriter(selected, false);
						writer.write("Monomer");
						writer.write(";");
						writer.write("mmol/gDW");
						writer.write(";");
						writer.write("Percentage (mol/mol)");
						writer.write("\r\n");
					for (String base : aaMolTextFieldMap.keySet()){
						writer.write(base);
						writer.write(";");
						writer.write(pContentmmol_gDW.get(base).toString());
						writer.write(";");
						writer.write(pContentm.get(base).toString());
						writer.write("\r\n");
					}
					System.out.println("Write success!");
					writer.close();
				} catch (IOException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
		        //System.out.print(pContentmmol_gDW.keySet()+ "-->" + pContentmmol_gDW.values()+ "-->" + pContentm.values());

			}      
		});
		
		// caixa de texto para inserir ficheiro com dados de express�o em .csv
		
		OpenGD.setBounds(251, 67, 72, 23);
		panel.add(OpenGD);
		
		JLabel lblExpressionCsvFile = new JLabel("Expression CSV file");
		lblExpressionCsvFile.setBounds(10, 71, 123, 14);
		panel.add(lblExpressionCsvFile);
		
		ExportP.setBounds(658, 98, 101, 23);
		panel.add(ExportP);
		
		// colocar as caixas de texto dos AA
		
		aaMolTextFieldMap = new HashMap<> ();
		for (int i = 0; i < aaArray.length / 2; i+=2) {
			String base = aaArray[i];
			String name = aaArray[i + 1];
			JTextField aaTextField = new JTextField();
			aaTextField.setColumns(10);
			
			JLabel aaLabel = new JLabel(name);
			
			panel_8.add(aaTextField);
			panel_8.add(aaLabel);
			
			aaMolTextFieldMap.put(base, aaTextField);
		}
		for (int i = aaArray.length / 2; i < aaArray.length; i+=2) {
			String base = aaArray[i];
			String name = aaArray[i + 1];
			JTextField aaTextField = new JTextField();
			aaTextField.setColumns(10);
			
			JLabel aaLabel = new JLabel(name);
			
			panel_7.add(aaTextField);
			panel_7.add(aaLabel);
			
			aaMolTextFieldMap.put(base, aaTextField);
		}

		

		panel_1 = new JPanel();
		tabbedPane.addTab("DNA", null, panel_1, null);
		panel_1.setLayout(null);

		JLabel lblDnaFastaFile = new JLabel("DNA FASTA file");
		lblDnaFastaFile.setBounds(10, 31, 85, 14);
		panel_1.add(lblDnaFastaFile);

		inputFD = new JTextField();
		inputFD.setColumns(10);
		inputFD.setBounds(128, 28, 113, 20);
		panel_1.add(inputFD);

		JButton OpenFD = new JButton("Open");
		OpenFD.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {

				FileDialog fd = new FileDialog(s, SWT.OPEN);
				fd.setText("Abrir");
				fd.setFilterPath("C:/");
				String[] filterExt = { "*.txt", "*.doc", ".rtf", "*.*" };
				fd.setFilterExtensions(filterExt);
				String selected = fd.open();
				inputFD.setText(selected);

			}
		});
		OpenFD.setBounds(251, 27, 72, 23);
		panel_1.add(OpenFD);

		JLabel lblDnaComposition = new JLabel("Cellular content DNA");
		lblDnaComposition.setBounds(518, 31, 119, 14);
		panel_1.add(lblDnaComposition);

		BM_dna = new JTextField();
		BM_dna.setColumns(10);
		BM_dna.setBounds(673, 28, 86, 20);
		panel_1.add(BM_dna);

		JLabel lblBiomassCompositionDntps = new JLabel(" Biomass Composition dNTPs (mmol/gDW)");
		lblBiomassCompositionDntps.setBounds(384, 144, 244, 14);
		panel_1.add(lblBiomassCompositionDntps);

		JPanel panel_5 = new JPanel();
		panel_5.setBounds(419, 178, 159, 108);
		panel_1.add(panel_5);

		A_DNA = new JTextField();
		A_DNA.setColumns(10);
		panel_5.add(A_DNA);

		JLabel lblDatp_1 = new JLabel("dATP");
		panel_5.add(lblDatp_1);

		T_DNA = new JTextField();
		T_DNA.setColumns(10);
		panel_5.add(T_DNA);

		JLabel lblDttp_1 = new JLabel("dTTP");
		panel_5.add(lblDttp_1);

		C_DNA = new JTextField();
		C_DNA.setColumns(10);
		panel_5.add(C_DNA);

		JLabel lblDctp = new JLabel("dCTP");
		panel_5.add(lblDctp);

		G_DNA = new JTextField();
		G_DNA.setColumns(10);
		panel_5.add(G_DNA);

		JLabel lblDgtp_1 = new JLabel("dGTP");
		panel_5.add(lblDgtp_1);

		JButton DetermineDNA = new JButton("Determine");
		DetermineDNA.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {

				try {

					DNAseq=(inputFD.getText());

					bm_dna=Double.parseDouble(BM_dna.getText());

					Object [] res=processaDNA();
					int tamanho=(Integer)res[1];
					Hashtable <NucleotideCompound, Integer> h1=(Hashtable<NucleotideCompound, Integer>) res[0];

					Hashtable <String, Float> MWDNA=new Hashtable <String, Float>();
					MWDNA.put("A",  (float) 491.2);
					MWDNA.put("C",  (float) 467.2);
					MWDNA.put("G",  (float) 507.2);
					MWDNA.put("T",  (float) 482.2);
					 

					Enumeration<NucleotideCompound> em=h1.keys();
					double soma=0;
					while (em.hasMoreElements()) {NucleotideCompound key1=em.nextElement();

					if(h1.get(key1)!=0 && key1.getShortName()!="S"){

						int count=h1.get(key1);
						double contentm_m=(float) count/(float)tamanho;
						dnaContentm.put(key1.getShortName(), contentm_m);
						double contentg=contentm_m*MWDNA.get(key1.getShortName());
						soma+=contentg;
						dnaContentg_m.put(key1.getShortName(), contentg);


						//	System.out.println(key1.getShortName()+"-->"+contentg_m.get(key1.getShortName()));

					}
					}

					double contentg_gA=dnaContentg_m.get("A")/soma;
					double contentg_gC=dnaContentg_m.get("C")/soma;
					double contentg_gG=dnaContentg_m.get("G")/soma;
					double contentg_gT=dnaContentg_m.get("T")/soma;
					


					double contentmmol_gDW_A=(float)(contentg_gA*bm_dna*1000)/MWDNA.get("A").shortValue();
					double contentmmol_gDW_C=(float)(contentg_gC*bm_dna*1000)/MWDNA.get("C").shortValue();
					double contentmmol_gDW_G=(float)(contentg_gG*bm_dna*1000)/MWDNA.get("G").shortValue();
					double contentmmol_gDW_T=(float)(contentg_gT*bm_dna*1000)/MWDNA.get("T").shortValue();
					
					dnaContentmmol_gDW.put("A", contentmmol_gDW_A);
					dnaContentmmol_gDW.put("C", contentmmol_gDW_C);
					dnaContentmmol_gDW.put("G", contentmmol_gDW_G);
					dnaContentmmol_gDW.put("T", contentmmol_gDW_T);
					
					dATP_mol.setText(String.format("%.5f", dnaContentm.get("A")));
					dTTP_mol.setText(String.format("%.5f", dnaContentm.get("T")));
					dCTP.setText(String.format("%.5f", dnaContentm.get("C")));
					dGTP.setText(String.format("%.5f", dnaContentm.get("G")));

					A_DNA.setText(String.format("%.5f", contentmmol_gDW_A));
					C_DNA.setText(String.format("%.5f", contentmmol_gDW_C));
					G_DNA.setText(String.format("%.5f", contentmmol_gDW_G));
					T_DNA.setText(String.format("%.5f", contentmmol_gDW_T));

				}



				catch (Exception e1) {e1.printStackTrace();}

			}
		});
		DetermineDNA.setBounds(658, 59, 101, 23);
		panel_1.add(DetermineDNA);

		JPanel panel_9 = new JPanel();
		panel_9.setBounds(177, 178, 159, 108);
		panel_1.add(panel_9);

		dATP_mol = new JTextField();
		dATP_mol.setColumns(10);
		panel_9.add(dATP_mol);

		JLabel lblDatp = new JLabel("dATP");
		panel_9.add(lblDatp);

		dTTP_mol = new JTextField();
		dTTP_mol.setColumns(10);
		panel_9.add(dTTP_mol);

		JLabel lblDttp = new JLabel("dTTP");
		panel_9.add(lblDttp);

		dCTP = new JTextField();
		dCTP.setColumns(10);
		panel_9.add(dCTP);

		JLabel lblCctp = new JLabel("dCTP");
		panel_9.add(lblCctp);

		dGTP = new JTextField();
		dGTP.setColumns(10);
		panel_9.add(dGTP);

		JLabel lblDgtp = new JLabel("dGTP");
		panel_9.add(lblDgtp);

		JLabel lblDntpContentmolmol = new JLabel("dNTP Content (mol/mol)");
		lblDntpContentmolmol.setBounds(201, 144, 148, 14);
		panel_1.add(lblDntpContentmolmol);
		
		JButton DExport = new JButton("Export");
		DExport.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				FileDialog gd = new FileDialog(s, SWT.SAVE);
				gd.setText("Guardar");
				gd.setFilterPath("C:/");
				String[] filterExt = { "*.csv"};
				gd.setFilterExtensions(filterExt);
		        String selected = gd.open();
		        gd.setFileName(selected);
		        fileName.setText(selected);
		        FileWriter writer;
		        try {
					writer = new FileWriter(selected, false);
						writer.write("Monomer");
						writer.write(";");
						writer.write("mmol/gDW");
						writer.write(";");
						writer.write("Percentage (mol/mol)");
						writer.write("\r\n");
					for (String base : dnaContentmmol_gDW.keySet()){
						writer.write(base);
						writer.write(";");
						writer.write(dnaContentmmol_gDW.get(base).toString());
						writer.write(";");
						writer.write(dnaContentm.get(base).toString());
						writer.write("\r\n");
					}
					System.out.println("Write success!");
					writer.close();
				} catch (IOException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
				
			}
		});
		DExport.setBounds(658, 93, 101, 23);
		panel_1.add(DExport);

		panel_2 = new JPanel();
		tabbedPane.addTab("RNA", null, panel_2, null);
		panel_2.setLayout(null);

		JLabel lblMrnaFastaFile = new JLabel("mRNA FASTA file");
		lblMrnaFastaFile.setBounds(10, 31, 101, 14);
		panel_2.add(lblMrnaFastaFile);

		inputFmRNA = new JTextField();
		inputFmRNA.setColumns(10);
		inputFmRNA.setBounds(128, 28, 113, 20);
		panel_2.add(inputFmRNA);

		JButton OpenmRNA = new JButton("Open");
		OpenmRNA.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {

				FileDialog fd = new FileDialog(s, SWT.OPEN);
				fd.setText("Abrir");
				fd.setFilterPath("C:/");
				String[] filterExt = { "*.txt", "*.doc", ".rtf", "*.*" };
				fd.setFilterExtensions(filterExt);
				String selected = fd.open();
				inputFmRNA.setText(selected);

			}
		});
		OpenmRNA.setBounds(251, 27, 72, 23);
		panel_2.add(OpenmRNA);

		JLabel lblTrnaFastaFile = new JLabel("tRNA FASTA file");
		lblTrnaFastaFile.setBounds(10, 59, 101, 14);
		panel_2.add(lblTrnaFastaFile);

		inputFtRNA = new JTextField();
		inputFtRNA.setColumns(10);
		inputFtRNA.setBounds(128, 56, 113, 20);
		panel_2.add(inputFtRNA);

		JButton OpentRNA = new JButton("Open");
		OpentRNA.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {

				FileDialog fd = new FileDialog(s, SWT.OPEN);
				fd.setText("Abrir");
				fd.setFilterPath("C:/");
				String[] filterExt = { "*.txt", "*.doc", ".rtf", "*.*" };
				fd.setFilterExtensions(filterExt);
				String selected = fd.open();
				inputFtRNA.setText(selected);

			}
		});
		OpentRNA.setBounds(251, 55, 72, 23);
		panel_2.add(OpentRNA);

		JLabel lblRrnaFastaFile = new JLabel("rRNA FASTA file");
		lblRrnaFastaFile.setBounds(10, 90, 101, 14);
		panel_2.add(lblRrnaFastaFile);

		inputFrRNA = new JTextField();
		inputFrRNA.setColumns(10);
		inputFrRNA.setBounds(128, 87, 113, 20);
		panel_2.add(inputFrRNA);

		JButton OpenrRNA = new JButton("Open");
		OpenrRNA.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {

				FileDialog fd = new FileDialog(s, SWT.OPEN);
				fd.setText("Abrir");
				fd.setFilterPath("C:/");
				String[] filterExt = { "*.txt", "*.doc", ".rtf", "*.*" };
				fd.setFilterExtensions(filterExt);
				String selected = fd.open();
				inputFrRNA.setText(selected);

			}
		});
		OpenrRNA.setBounds(251, 86, 72, 23);
		panel_2.add(OpenrRNA);

		JLabel lblRnaComposition = new JLabel("Cellular content RNA");
		lblRnaComposition.setBounds(518, 31, 145, 14);
		panel_2.add(lblRnaComposition);

		BM_rna = new JTextField();
		BM_rna.setColumns(10);
		BM_rna.setBounds(673, 28, 86, 20);
		panel_2.add(BM_rna);

		JLabel lblBiomassCompositon = new JLabel("Biomass Compositon NTP (mmol/gDW) ");
		lblBiomassCompositon.setBounds(400, 144, 275, 14);
		panel_2.add(lblBiomassCompositon);

		JButton DetermineRNA = new JButton("Determine");
		DetermineRNA.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {

				try {

					mRNAseq=(inputFmRNA.getText());
					perc_mRNA=Double.parseDouble(mRNA_per.getText());
					
					tRNAseq=(inputFtRNA.getText());
					perc_tRNA=Double.parseDouble(tRNA_per.getText());
					
					rRNAseq=(inputFrRNA.getText());
					perc_rRNA=Double.parseDouble(rRNA_per.getText());

					Object [] resm=processamRNA();
					int tamanhomRNA=(Integer)resm[1];
					
					Object [] rest=processatRNA();
					int tamanhotRNA=(Integer)rest[1];
					
					Object [] resr=processarRNA();
					int tamanhorRNA=(Integer)resr[1];

					Hashtable <NucleotideCompound, Integer> h1m = (Hashtable <NucleotideCompound, Integer>)resm[0];
					Hashtable <NucleotideCompound, Integer> h1t = (Hashtable<NucleotideCompound, Integer>) rest[0];
					Hashtable <NucleotideCompound, Integer> h1r=(Hashtable<NucleotideCompound, Integer>) resr[0];

					Enumeration<NucleotideCompound> emm=h1m.keys();
					while (emm.hasMoreElements()) {NucleotideCompound key1=emm.nextElement();

					if(h1m.get(key1)!=0){

						int countmRNA=h1m.get(key1);
						double contentm_m1=(float) (perc_mRNA*countmRNA/tamanhomRNA);
						rnaContentm_mm.put(key1.getShortName(), contentm_m1);
						
						int counttRNA=h1t.get(key1);
						double contentm_m2=(float) (perc_tRNA*counttRNA/tamanhotRNA);
						rnaContentm_mt.put(key1.getShortName(), contentm_m2);
						
						int countrRNA=h1r.get(key1);			
						double contentm_m3=(float) (perc_rRNA*countrRNA/tamanhorRNA);
						rnaContentm_mr.put(key1.getShortName(), contentm_m3);
						
						rnaContentm_m.put(key1.getShortName(), rnaContentm_mm.get(key1.getShortName())+ rnaContentm_mt.get(key1.getShortName())+rnaContentm_mr.get(key1.getShortName()));
						
						
				}
					}

					rnaContentm_m.put("U", rnaContentm_m.get("T"));
					
					ATP_mol.setText(String.format("%.5f", rnaContentm_m.get("A")));
					CTP_mol.setText(String.format("%.5f", rnaContentm_m.get("C")));
					GTP_mol.setText(String.format("%.5f", rnaContentm_m.get("G")));
					UTP_mol.setText(String.format("%.5f", rnaContentm_m.get("T")));
					/*
					Hashtable <String, Float> MWRNA=new Hashtable <String, Float>();
					MWRNA.put("A",  (float) 328.201);
					MWRNA.put("C",  (float) 304.175);
					MWRNA.put("G",  (float) 344.2);
					MWRNA.put("U",  (float) 305.201);
*/
					
					Hashtable <String, Float> MWRNA=new Hashtable <String, Float>();
					MWRNA.put("A",  (float) 503.2);
					MWRNA.put("C",  (float) 479.1);
					MWRNA.put("G",  (float) 519.1);
					MWRNA.put("U",  (float) 480.1);
					 

					double contentg_mA=rnaContentm_m.get("A")*MWRNA.get("A");
					double contentg_mU=rnaContentm_m.get("T")*MWRNA.get("U");
					double contentg_mC=rnaContentm_m.get("C")*MWRNA.get("C");
					double contentg_mG=rnaContentm_m.get("G")*MWRNA.get("G");

					double soma=contentg_mA+contentg_mU+contentg_mC+contentg_mG;

					double contentg_gA=contentg_mA/soma;
					double contentg_gC=contentg_mC/soma;
					double contentg_gG=contentg_mG/soma;
					double contentg_gU=contentg_mU/soma;

					bm_rna=Double.parseDouble(BM_rna.getText());

					double rnaContentmmol_gDW_A=(double)(contentg_gA*bm_rna*1000)/MWRNA.get("A").shortValue();
					double rnaContentmmol_gDW_C=(double)(contentg_gC*bm_rna*1000)/MWRNA.get("C").shortValue();
					double rnaContentmmol_gDW_G=(double)(contentg_gG*bm_rna*1000)/MWRNA.get("G").shortValue();
					double rnaContentmmol_gDW_U=(double)(contentg_gU*bm_rna*1000)/MWRNA.get("U").shortValue();

					rnaContentmmol_gDW.put("A", rnaContentmmol_gDW_A);
					rnaContentmmol_gDW.put("C", rnaContentmmol_gDW_C);
					rnaContentmmol_gDW.put("G", rnaContentmmol_gDW_G);
					rnaContentmmol_gDW.put("U", rnaContentmmol_gDW_U);
					
					RNA_A.setText(String.format("%.5f", rnaContentmmol_gDW_A));
					RNA_C.setText(String.format("%.5f", rnaContentmmol_gDW_C));
					RNA_G.setText(String.format("%.5f", rnaContentmmol_gDW_G));
					RNA_U.setText(String.format("%.5f",rnaContentmmol_gDW_U));

				}

				catch (Exception e1) {e1.printStackTrace();}

			}
		});
		DetermineRNA.setBounds(658, 59, 101, 23);
		panel_2.add(DetermineRNA);

		JPanel panel_6 = new JPanel();
		panel_6.setBounds(419, 178, 157, 108);
		panel_2.add(panel_6);

		RNA_A = new JTextField();
		RNA_A.setColumns(10);
		panel_6.add(RNA_A);

		JLabel lblAtp_1 = new JLabel("ATP");
		panel_6.add(lblAtp_1);

		RNA_U = new JTextField();
		RNA_U.setColumns(10);
		panel_6.add(RNA_U);

		JLabel lblUtp_1 = new JLabel("UTP");
		panel_6.add(lblUtp_1);

		RNA_C = new JTextField();
		RNA_C.setColumns(10);
		panel_6.add(RNA_C);

		JLabel lblCtp_1 = new JLabel("CTP");
		panel_6.add(lblCtp_1);

		RNA_G = new JTextField();
		RNA_G.setColumns(10);
		panel_6.add(RNA_G);

		JLabel lblGtp_1 = new JLabel("GTP");
		panel_6.add(lblGtp_1);

		JLabel lblNtpContent = new JLabel("NTP Content (mol/mol)");
		lblNtpContent.setBounds(201, 144, 192, 14);
		panel_2.add(lblNtpContent);

		JPanel panel_10 = new JPanel();
		panel_10.setBounds(177, 178, 157, 108);
		panel_2.add(panel_10);

		ATP_mol = new JTextField();
		ATP_mol.setColumns(10);
		panel_10.add(ATP_mol);

		JLabel lblAtp = new JLabel("ATP");
		panel_10.add(lblAtp);

		UTP_mol = new JTextField();
		UTP_mol.setColumns(10);
		panel_10.add(UTP_mol);

		JLabel lblUtp = new JLabel("UTP");
		panel_10.add(lblUtp);

		CTP_mol = new JTextField();
		CTP_mol.setColumns(10);
		panel_10.add(CTP_mol);

		JLabel lblCtp = new JLabel("CTP");
		panel_10.add(lblCtp);

		GTP_mol = new JTextField();
		GTP_mol.setColumns(10);
		panel_10.add(GTP_mol);

		JLabel lblGtp = new JLabel("GTP");
		panel_10.add(lblGtp);

		JLabel lblMrnaPercentage = new JLabel("mRNA percentage");
		lblMrnaPercentage.setBounds(341, 31, 113, 14);
		panel_2.add(lblMrnaPercentage);

		mRNA_per = new JTextField();
		mRNA_per.setColumns(10);
		mRNA_per.setBounds(464, 28, 44, 20);
		panel_2.add(mRNA_per);

		JLabel lblTrnaPercentage = new JLabel("tRNA percentage");
		lblTrnaPercentage.setBounds(341, 59, 113, 14);
		panel_2.add(lblTrnaPercentage);

		JLabel lblRrnaPercentage = new JLabel("rRNA percentage");
		lblRrnaPercentage.setBounds(341, 90, 113, 14);
		panel_2.add(lblRrnaPercentage);

		tRNA_per = new JTextField();
		tRNA_per.setColumns(10);
		tRNA_per.setBounds(464, 56, 44, 20);
		panel_2.add(tRNA_per);

		rRNA_per = new JTextField();
		rRNA_per.setColumns(10);
		rRNA_per.setBounds(464, 87, 44, 20);
		panel_2.add(rRNA_per);
		
		JButton RExport = new JButton("Export");
		RExport.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				FileDialog gd = new FileDialog(s, SWT.SAVE);
				gd.setText("Guardar");
				gd.setFilterPath("C:/");
				String[] filterExt = { "*.csv"};
				gd.setFilterExtensions(filterExt);
		        String selected = gd.open();
		        gd.setFileName(selected);
		        fileName.setText(selected);
		        FileWriter writer;
		        try {
					writer = new FileWriter(selected, false);
						writer.write("Monomer");
						writer.write(";");
						writer.write("mmol/gDW");
						writer.write(";");
						writer.write("Percentage (mol/mol)");
						writer.write("\r\n");
					for (String base : rnaContentmmol_gDW.keySet()){
						writer.write(base);
						writer.write(";");
						writer.write(rnaContentmmol_gDW.get(base).toString());
						writer.write(";");
						writer.write(rnaContentm_m.get(base).toString());
						writer.write("\r\n");
					}
					System.out.println("Write success!");
					writer.close();
				} catch (IOException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
				
			}
		});
		RExport.setBounds(658, 90, 101, 23);
		panel_2.add(RExport);



	}
	
	private Object [] processaProteinas ()throws Exception {
		File file= new File (AAseq); // ir buscar caminho para o ficheiro � caixa DNAseq
		int tamanho=0;//n total de amino�cidos
		Object [] ret=new Object[2];


		Hashtable <AminoAcidCompound, Integer> h=new Hashtable <AminoAcidCompound, Integer>(); //Hash com o AA e o numero de cada amino�cido
		Map<String, ProteinSequence> sequencias=FastaReaderHelper.readFastaProteinSequence(file); //ler ficheiro Fasta

		Iterator<String> it=sequencias.keySet().iterator(); //iterador pelas sequencias 
		while (it.hasNext()) {

			String key=it.next();
			ProteinSequence protseq=sequencias.get(key);
			tamanho+=protseq.getLength();

			List<AminoAcidCompound> l=protseq.getCompoundSet().getAllCompounds();
			for (int i=0;i<l.size();i++) {
				int count=0;
				AminoAcidCompound C=l.get(i);
				if(h.containsKey(C)) {
					count=protseq.countCompounds(l.get(i))+h.get(C);}

				else {
					count=protseq.countCompounds(l.get(i));
				}
				h.put(C, count);

			}
			//					System.out.print(protID);

		}
		tamanho=tamanho+sequencias.size();

		ret[0]=h;
		ret[1]=new Integer (tamanho);
		return ret;

	}

	private Object [] processaDNA ()throws Exception {
		File file= new File (DNAseq);
		int tamanho=0;

		Object [] ret=new Object[2];
		Hashtable <NucleotideCompound, Integer> h=new Hashtable <NucleotideCompound, Integer>(); //Hash com o AA e o numero de cada amino�cido

		LinkedHashMap <String, DNASequence> sequenciasDNA=FastaReaderHelper.readFastaDNASequence(file); //ler ficheiro Fasta
		Iterator<String> it=sequenciasDNA.keySet().iterator(); //iterador pelas sequencias
		while (it.hasNext()) {
			String key=it.next();
			DNASequence dnaseq=sequenciasDNA.get(key);
			tamanho+=dnaseq.getLength();

			List<NucleotideCompound> l=dnaseq.getCompoundSet().getAllCompounds();
			for (int i=0;i<l.size();i++) {
				int count=0;
				NucleotideCompound C=l.get(i);
				if(h.containsKey(C)) {
					count=dnaseq.countCompounds(l.get(i))+h.get(C);}
				else {
					count=dnaseq.countCompounds(l.get(i));
				}
				h.put(C, count);
			}

		}

		ret[0]=h;
		ret[1]=new Integer (tamanho);

		return ret;
	}

	private Object [] processamRNA ()throws Exception {
		File filemRNA= new File (mRNAseq);

		int tamanhom=0;
		int countm=0;

		Object [] retm=new Object[2];

		Hashtable <NucleotideCompound, Integer> h=new Hashtable <NucleotideCompound, Integer>(); //Hash com o AA e o numero de cada amino�cido
		LinkedHashMap <String, DNASequence> sequenciasmRNA=FastaReaderHelper.readFastaDNASequence(filemRNA);
		Iterator<String> itm=sequenciasmRNA.keySet().iterator(); //iterador pelas sequencias

		while (itm.hasNext()) {
			String key=itm.next();
			DNASequence rnaseq=sequenciasmRNA.get(key);
			tamanhom+=rnaseq.getLength();

			List<NucleotideCompound> lm=rnaseq.getCompoundSet().getAllCompounds();
			for (int i=0;i<lm.size();i++) {
				countm=0;
				NucleotideCompound C=lm.get(i);
				if(h.containsKey(C)) {
					countm=rnaseq.countCompounds(lm.get(i))+h.get(C);}
				else {
					countm=rnaseq.countCompounds(lm.get(i));
				}
				h.put(C, countm);
			}

		}

		retm[0]=h;
		retm[1]=new Integer (tamanhom);

		return retm;
	}

	private Object [] processatRNA ()throws Exception {
		File filetRNA= new File (tRNAseq);

		int tamanhot=0;

		Object [] rett=new Object[2];

		Hashtable <NucleotideCompound, Integer> h=new Hashtable <NucleotideCompound, Integer>(); //Hash com o AA e o numero de cada amino�cido
		LinkedHashMap <String, DNASequence> sequenciastRNA=FastaReaderHelper.readFastaDNASequence(filetRNA);
		Iterator<String> itt=sequenciastRNA.keySet().iterator(); //iterador pelas sequencias

		while (itt.hasNext()) {
			String key=itt.next();
			DNASequence rnaseq=sequenciastRNA.get(key);
			tamanhot+=rnaseq.getLength();

			List<NucleotideCompound> lt=rnaseq.getCompoundSet().getAllCompounds();
			for (int i=0;i<lt.size();i++) {
				int countt=0;
				NucleotideCompound C=lt.get(i);
				if(h.containsKey(C)) {
					countt=rnaseq.countCompounds(lt.get(i))+h.get(C);}
				else {
					countt=rnaseq.countCompounds(lt.get(i));
				}
				h.put(C, countt);
			}

		}

		rett[0]=h;
		rett[1]=new Integer (tamanhot);
		System.out.print (tamanhot);
		return rett;
	}

	private Object [] processarRNA ()throws Exception {
		File filerRNA= new File (rRNAseq);

		int tamanhor=0;

		Object [] retr=new Object[2];

		Hashtable <NucleotideCompound, Integer> h=new Hashtable <NucleotideCompound, Integer>(); //Hash com o AA e o numero de cada amino�cido
		LinkedHashMap <String, DNASequence> sequenciasrRNA=FastaReaderHelper.readFastaDNASequence(filerRNA);
		Iterator<String> itr=sequenciasrRNA.keySet().iterator(); //iterador pelas sequencias

		while (itr.hasNext()) {
			String key=itr.next();
			DNASequence rnaseq=sequenciasrRNA.get(key);
			tamanhor+=rnaseq.getLength();

			List<NucleotideCompound> lr=rnaseq.getCompoundSet().getAllCompounds();
			for (int i=0;i<lr.size();i++) {
				int countr=0;
				NucleotideCompound C=lr.get(i);
				if(h.containsKey(C)) {
					countr=rnaseq.countCompounds(lr.get(i))+h.get(C);}
				else {
					countr=rnaseq.countCompounds(lr.get(i));
				}
				h.put(C, countr);
			}

		}

		retr[0]=h;
		retr[1]=new Integer (tamanhor);

		return retr;
	}
}
