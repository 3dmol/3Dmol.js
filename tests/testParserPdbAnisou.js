viewer = $3Dmol.createViewer("#gldiv");

      var data = $("#pdbData").val(); 
      //console.log(data);

      viewer.setBackgroundColor(0xffffff);
      viewer.addModel(data, "pdb");
      viewer.zoomTo();
      viewer.render();

      var model = viewer.getModel();
      var atoms = model.selectedAtoms();
      var html = "";
      for (var i = 0; i < atoms.length; i++) {
        var a = atoms[i];
        if ("uMat" in a) {
          html += a.serial + " " + JSON.stringify(a["uMat"]) + "<br/>";
        }
      }
      $("#loadAnisouResultDiv").html(html);
      console.log(atoms[0]);

/* @data


HEADER    ISOMERASE                               25-FEB-14   4POC              
TITLE     STRUCTURE OF TRIOSEPHOSPHATE ISOMERASE WILD TYPE HUMAN ENZYME.        
COMPND    MOL_ID: 1;                                                            
COMPND   2 MOLECULE: TRIOSEPHOSPHATE ISOMERASE;                                 
COMPND   3 CHAIN: A, B;                                                         
COMPND   4 SYNONYM: TIM, TRIOSE-PHOSPHATE ISOMERASE;                            
COMPND   5 EC: 5.3.1.1;                                                         
COMPND   6 ENGINEERED: YES                                                      
SOURCE    MOL_ID: 1;                                                            
SOURCE   2 ORGANISM_SCIENTIFIC: HOMO SAPIENS;                                   
SOURCE   3 ORGANISM_COMMON: HUMAN;                                              
SOURCE   4 ORGANISM_TAXID: 9606;                                                
SOURCE   5 GENE: TPI, TPI1;                                                     
SOURCE   6 EXPRESSION_SYSTEM: ESCHERICHIA COLI;                                 
SOURCE   7 EXPRESSION_SYSTEM_TAXID: 469008;                                     
SOURCE   8 EXPRESSION_SYSTEM_STRAIN: BL21(DE3) CODON+RILP;                      
SOURCE   9 EXPRESSION_SYSTEM_VECTOR_TYPE: PLASMID;                              
SOURCE  10 EXPRESSION_SYSTEM_PLASMID: PLC3                                      
KEYWDS    TIM ALPHA/BETA BARREL, TIM BARREL, ISOMERASE, GLYCOLYTIC              
EXPDTA    X-RAY DIFFRACTION                                                     
AUTHOR    C.G.AMRICH,A.A.ASLAM,A.HEROUX,A.P.VANDEMARK                           
REVDAT   1   14-JAN-15 4POC    0                                                
JRNL        AUTH   B.P.ROLAND,C.G.AMRICH,C.J.KAMMERER,K.A.STUCHUL,S.B.LARSEN,   
JRNL        AUTH 2 S.RODE,A.A.ASLAM,A.HEROUX,R.WETZEL,A.P.VANDEMARK,            
JRNL        AUTH 3 M.J.PALLADINO                                                
JRNL        TITL   TRIOSEPHOSPHATE ISOMERASE I170V ALTERS CATALYTIC SITE,       
JRNL        TITL 2 ENHANCES STABILITY AND INDUCES PATHOLOGY IN A DROSOPHILA     
JRNL        TITL 3 MODEL OF TPI DEFICIENCY.                                     
JRNL        REF    BIOCHIM.BIOPHYS.ACTA          V.1852    61 2015              
JRNL        REFN                   ISSN 0006-3002                               
JRNL        PMID   25463631                                                     
JRNL        DOI    10.1016/J.BBADIS.2014.10.010                                 
REMARK   2                                                                      
REMARK   2 RESOLUTION.    1.60 ANGSTROMS.                                       
REMARK   3                                                                      
REMARK   3 REFINEMENT.                                                          
REMARK   3   PROGRAM     : PHENIX (PHENIX.REFINE: 1.8.2_1309)                   
REMARK   3   AUTHORS     : PAUL ADAMS,PAVEL AFONINE,VINCENT CHEN,IAN            
REMARK   3               : DAVIS,KRESHNA GOPAL,RALF GROSSE-                     
REMARK   3               : KUNSTLEVE,LI-WEI HUNG,ROBERT IMMORMINO,              
REMARK   3               : TOM IOERGER,AIRLIE MCCOY,ERIK MCKEE,NIGEL            
REMARK   3               : MORIARTY,REETAL PAI,RANDY READ,JANE                  
REMARK   3               : RICHARDSON,DAVID RICHARDSON,TOD ROMO,JIM             
REMARK   3               : SACCHETTINI,NICHOLAS SAUTER,JACOB SMITH,             
REMARK   3               : LAURENT STORONI,TOM TERWILLIGER,PETER                
REMARK   3               : ZWART                                                
REMARK   3                                                                      
REMARK   3    REFINEMENT TARGET : ML                                            
REMARK   3                                                                      
REMARK   3  DATA USED IN REFINEMENT.                                            
REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : 1.60                           
REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : 46.56                          
REMARK   3   MIN(FOBS/SIGMA_FOBS)              : 0.000                          
REMARK   3   COMPLETENESS FOR RANGE        (%) : 93.7                           
REMARK   3   NUMBER OF REFLECTIONS             : 52340                          
REMARK   3                                                                      
REMARK   3  FIT TO DATA USED IN REFINEMENT.                                     
REMARK   3   R VALUE     (WORKING + TEST SET) : 0.155                           
REMARK   3   R VALUE            (WORKING SET) : 0.153                           
REMARK   3   FREE R VALUE                     : 0.187                           
REMARK   3   FREE R VALUE TEST SET SIZE   (%) : 5.150                           
REMARK   3   FREE R VALUE TEST SET COUNT      : 2695                            
REMARK   3                                                                      
REMARK   3  FIT TO DATA USED IN REFINEMENT (IN BINS).                           
REMARK   3   BIN  RESOLUTION RANGE  COMPL.    NWORK NFREE   RWORK  RFREE        
REMARK   3     1 46.5840 -  4.2713    0.99     2859   153  0.1558 0.1716        
REMARK   3     2  4.2713 -  3.3906    0.99     2793   145  0.1393 0.1563        
REMARK   3     3  3.3906 -  2.9620    0.99     2783   154  0.1523 0.1960        
REMARK   3     4  2.9620 -  2.6913    1.00     2816   150  0.1524 0.1852        
REMARK   3     5  2.6913 -  2.4984    1.00     2771   163  0.1530 0.2025        
REMARK   3     6  2.4984 -  2.3511    1.00     2784   145  0.1501 0.2088        
REMARK   3     7  2.3511 -  2.2333    1.00     2826   121  0.1433 0.1765        
REMARK   3     8  2.2333 -  2.1361    1.00     2763   164  0.1438 0.1789        
REMARK   3     9  2.1361 -  2.0539    1.00     2766   167  0.1609 0.2152        
REMARK   3    10  2.0539 -  1.9830    1.00     2773   165  0.1594 0.1842        
REMARK   3    11  1.9830 -  1.9210    1.00     2766   137  0.1560 0.2126        
REMARK   3    12  1.9210 -  1.8661    1.00     2797   146  0.1530 0.2058        
REMARK   3    13  1.8661 -  1.8170    1.00     2788   142  0.1592 0.2306        
REMARK   3    14  1.8170 -  1.7726    1.00     2782   137  0.1706 0.2331        
REMARK   3    15  1.7726 -  1.7323    1.00     2766   156  0.1737 0.2026        
REMARK   3    16  1.7323 -  1.6955    0.91     2479   144  0.1808 0.2305        
REMARK   3    17  1.6955 -  1.6615    0.76     2097   134  0.1960 0.2225        
REMARK   3    18  1.6615 -  1.6302    0.64     1747   104  0.2237 0.2766        
REMARK   3    19  1.6302 -  1.6011    0.53     1489    68  0.2425 0.3406        
REMARK   3                                                                      
REMARK   3  BULK SOLVENT MODELLING.                                             
REMARK   3   METHOD USED        : FLAT BULK SOLVENT MODEL                       
REMARK   3   SOLVENT RADIUS     : 1.11                                          
REMARK   3   SHRINKAGE RADIUS   : 0.90                                          
REMARK   3   K_SOL              : NULL                                          
REMARK   3   B_SOL              : NULL                                          
REMARK   3                                                                      
REMARK   3  ERROR ESTIMATES.                                                    
REMARK   3   COORDINATE ERROR (MAXIMUM-LIKELIHOOD BASED)     : 0.190            
REMARK   3   PHASE ERROR (DEGREES, MAXIMUM-LIKELIHOOD BASED) : 20.720           
REMARK   3                                                                      
REMARK   3  B VALUES.                                                           
REMARK   3   FROM WILSON PLOT           (A**2) : 25.84                          
REMARK   3   MEAN B VALUE      (OVERALL, A**2) : 36.13                          
REMARK   3   OVERALL ANISOTROPIC B VALUE.                                       
REMARK   3    B11 (A**2) : NULL                                                 
REMARK   3    B22 (A**2) : NULL                                                 
REMARK   3    B33 (A**2) : NULL                                                 
REMARK   3    B12 (A**2) : NULL                                                 
REMARK   3    B13 (A**2) : NULL                                                 
REMARK   3    B23 (A**2) : NULL                                                 
REMARK   3                                                                      
REMARK   3  TWINNING INFORMATION.                                               
REMARK   3   FRACTION: NULL                                                     
REMARK   3   OPERATOR: NULL                                                     
REMARK   3                                                                      
REMARK   3  DEVIATIONS FROM IDEAL VALUES.                                       
REMARK   3                 RMSD          COUNT                                  
REMARK   3   BOND      :  0.005           3792                                  
REMARK   3   ANGLE     :  0.962           5137                                  
REMARK   3   CHIRALITY :  0.052            581                                  
REMARK   3   PLANARITY :  0.004            665                                  
REMARK   3   DIHEDRAL  : 15.111           1374                                  
REMARK   3                                                                      
REMARK   3  TLS DETAILS                                                         
REMARK   3   NUMBER OF TLS GROUPS  : NULL                                       
REMARK   3                                                                      
REMARK   3  NCS DETAILS                                                         
REMARK   3   NUMBER OF NCS GROUPS : NULL                                        
REMARK   3                                                                      
REMARK   3  OTHER REFINEMENT REMARKS: NULL                                      
REMARK   4                                                                      
REMARK   4 4POC COMPLIES WITH FORMAT V. 3.30, 13-JUL-11                         
REMARK 100                                                                      
REMARK 100 THIS ENTRY HAS BEEN PROCESSED BY RCSB ON 26-FEB-14.                  
REMARK 100 THE RCSB ID CODE IS RCSB085018.                                      
REMARK 200                                                                      
REMARK 200 EXPERIMENTAL DETAILS                                                 
REMARK 200  EXPERIMENT TYPE                : X-RAY DIFFRACTION                  
REMARK 200  DATE OF DATA COLLECTION        : 23-MAR-12                          
REMARK 200  TEMPERATURE           (KELVIN) : 100                                
REMARK 200  PH                             : 7.5                                
REMARK 200  NUMBER OF CRYSTALS USED        : 1                                  
REMARK 200                                                                      
REMARK 200  SYNCHROTRON              (Y/N) : Y                                  
REMARK 200  RADIATION SOURCE               : NSLS                               
REMARK 200  BEAMLINE                       : X25                                
REMARK 200  X-RAY GENERATOR MODEL          : NULL                               
REMARK 200  MONOCHROMATIC OR LAUE    (M/L) : M                                  
REMARK 200  WAVELENGTH OR RANGE        (A) : 1.1                                
REMARK 200  MONOCHROMATOR                  : SI-111                             
REMARK 200  OPTICS                         : NULL                               
REMARK 200                                                                      
REMARK 200  DETECTOR TYPE                  : CCD                                
REMARK 200  DETECTOR MANUFACTURER          : ADSC QUANTUM 315                   
REMARK 200  INTENSITY-INTEGRATION SOFTWARE : DENZO                              
REMARK 200  DATA SCALING SOFTWARE          : SCALEPACK                          
REMARK 200                                                                      
REMARK 200  NUMBER OF UNIQUE REFLECTIONS   : 52380                              
REMARK 200  RESOLUTION RANGE HIGH      (A) : 1.600                              
REMARK 200  RESOLUTION RANGE LOW       (A) : 50.000                             
REMARK 200  REJECTION CRITERIA  (SIGMA(I)) : 2.000                              
REMARK 200                                                                      
REMARK 200 OVERALL.                                                             
REMARK 200  COMPLETENESS FOR RANGE     (%) : 93.8                               
REMARK 200  DATA REDUNDANCY                : 5.600                              
REMARK 200  R MERGE                    (I) : 0.07700                            
REMARK 200  R SYM                      (I) : NULL                               
REMARK 200  <I/SIGMA(I)> FOR THE DATA SET  : 11.3000                            
REMARK 200                                                                      
REMARK 200 IN THE HIGHEST RESOLUTION SHELL.                                     
REMARK 200  HIGHEST RESOLUTION SHELL, RANGE HIGH (A) : 1.60                     
REMARK 200  HIGHEST RESOLUTION SHELL, RANGE LOW  (A) : 1.63                     
REMARK 200  COMPLETENESS FOR SHELL     (%) : 54.2                               
REMARK 200  DATA REDUNDANCY IN SHELL       : 3.00                               
REMARK 200  R MERGE FOR SHELL          (I) : 0.44200                            
REMARK 200  R SYM FOR SHELL            (I) : NULL                               
REMARK 200  <I/SIGMA(I)> FOR SHELL         : NULL                               
REMARK 200                                                                      
REMARK 200 DIFFRACTION PROTOCOL: SINGLE WAVELENGTH                              
REMARK 200 METHOD USED TO DETERMINE THE STRUCTURE: MOLECULAR REPLACEMENT        
REMARK 200 SOFTWARE USED: PHASER                                                
REMARK 200 STARTING MODEL: PDB ENTRY 2JK2                                       
REMARK 200                                                                      
REMARK 200 REMARK: NULL                                                         
REMARK 280                                                                      
REMARK 280 CRYSTAL                                                              
REMARK 280 SOLVENT CONTENT, VS   (%): 37.33                                     
REMARK 280 MATTHEWS COEFFICIENT, VM (ANGSTROMS**3/DA): 1.96                     
REMARK 280                                                                      
REMARK 280 CRYSTALLIZATION CONDITIONS: 35% PEG 2000 MME, 0.05 KBR, PH 7.5,      
REMARK 280  VAPOR DIFFUSION, SITTING DROP, TEMPERATURE 277K                     
REMARK 290                                                                      
REMARK 290 CRYSTALLOGRAPHIC SYMMETRY                                            
REMARK 290 SYMMETRY OPERATORS FOR SPACE GROUP: P 1 21 1                         
REMARK 290                                                                      
REMARK 290      SYMOP   SYMMETRY                                                
REMARK 290     NNNMMM   OPERATOR                                                
REMARK 290       1555   X,Y,Z                                                   
REMARK 290       2555   -X,Y+1/2,-Z                                             
REMARK 290                                                                      
REMARK 290     WHERE NNN -> OPERATOR NUMBER                                     
REMARK 290           MMM -> TRANSLATION VECTOR                                  
REMARK 290                                                                      
REMARK 290 CRYSTALLOGRAPHIC SYMMETRY TRANSFORMATIONS                            
REMARK 290 THE FOLLOWING TRANSFORMATIONS OPERATE ON THE ATOM/HETATM             
REMARK 290 RECORDS IN THIS ENTRY TO PRODUCE CRYSTALLOGRAPHICALLY                
REMARK 290 RELATED MOLECULES.                                                   
REMARK 290   SMTRY1   1  1.000000  0.000000  0.000000        0.00000            
REMARK 290   SMTRY2   1  0.000000  1.000000  0.000000        0.00000            
REMARK 290   SMTRY3   1  0.000000  0.000000  1.000000        0.00000            
REMARK 290   SMTRY1   2 -1.000000  0.000000  0.000000        0.00000            
REMARK 290   SMTRY2   2  0.000000  1.000000  0.000000       24.42550            
REMARK 290   SMTRY3   2  0.000000  0.000000 -1.000000        0.00000            
REMARK 290                                                                      
REMARK 290 REMARK: NULL                                                         
REMARK 300                                                                      
REMARK 300 BIOMOLECULE: 1                                                       
REMARK 300 SEE REMARK 350 FOR THE AUTHOR PROVIDED AND/OR PROGRAM                
REMARK 300 GENERATED ASSEMBLY INFORMATION FOR THE STRUCTURE IN                  
REMARK 300 THIS ENTRY. THE REMARK MAY ALSO PROVIDE INFORMATION ON               
REMARK 300 BURIED SURFACE AREA.                                                 
REMARK 350                                                                      
REMARK 350 COORDINATES FOR A COMPLETE MULTIMER REPRESENTING THE KNOWN           
REMARK 350 BIOLOGICALLY SIGNIFICANT OLIGOMERIZATION STATE OF THE                
REMARK 350 MOLECULE CAN BE GENERATED BY APPLYING BIOMT TRANSFORMATIONS          
REMARK 350 GIVEN BELOW.  BOTH NON-CRYSTALLOGRAPHIC AND                          
REMARK 350 CRYSTALLOGRAPHIC OPERATIONS ARE GIVEN.                               
REMARK 350                                                                      
REMARK 350 BIOMOLECULE: 1                                                       
REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: DIMERIC                           
REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE: DIMERIC                    
REMARK 350 SOFTWARE USED: PISA                                                  
REMARK 350 TOTAL BURIED SURFACE AREA: 4530 ANGSTROM**2                          
REMARK 350 SURFACE AREA OF THE COMPLEX: 18810 ANGSTROM**2                       
REMARK 350 CHANGE IN SOLVENT FREE ENERGY: -38.0 KCAL/MOL                        
REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, B                                  
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000            
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000            
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000            
REMARK 465                                                                      
REMARK 465 MISSING RESIDUES                                                     
REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE                       
REMARK 465 EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN               
REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)                
REMARK 465                                                                      
REMARK 465   M RES C SSSEQI                                                     
REMARK 465     Gly A    -5                                                      
REMARK 465     ASP A    -4                                                      
REMARK 465     Ile A    -3                                                      
REMARK 465     Thr A    -2                                                      
REMARK 465     HIS A    -1                                                      
REMARK 465     MET A     0                                                      
REMARK 465     Ala A     1                                                      
REMARK 465     PRO A     2                                                      
REMARK 465     Gly B    -5                                                      
REMARK 465     ASP B    -4                                                      
REMARK 465     Ile B    -3                                                      
REMARK 465     Thr B    -2                                                      
REMARK 465     HIS B    -1                                                      
REMARK 465     MET B     0                                                      
REMARK 465     Ala B     1                                                      
REMARK 500                                                                      
REMARK 500 GEOMETRY AND STEREOCHEMISTRY                                         
REMARK 500 SUBTOPIC: TORSION ANGLES                                             
REMARK 500                                                                      
REMARK 500 TORSION ANGLES OUTSIDE THE EXPECTED RAMACHANDRAN REGIONS:            
REMARK 500 (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN IDENTIFIER;               
REMARK 500 SSEQ=SEQUENCE NUMBER; I=INSERTION CODE).                             
REMARK 500                                                                      
REMARK 500 STANDARD TABLE:                                                      
REMARK 500 FORMAT:(10X,I3,1X,A3,1X,A1,I4,A1,4X,F7.2,3X,F7.2)                    
REMARK 500                                                                      
REMARK 500 EXPECTED VALUES: GJ KLEYWEGT AND TA JONES (1996). PHI/PSI-           
REMARK 500 CHOLOGY: RAMACHANDRAN REVISITED. STRUCTURE 4, 1395 - 1400            
REMARK 500                                                                      
REMARK 500  M RES CSSEQI        PSI       PHI                                   
REMARK 500    LYS A  13     -149.26     54.18                                   
REMARK 500    VAL A 196      -79.68   -121.56                                   
REMARK 500    LYS A 247       69.58   -101.29                                   
REMARK 500    LYS B  13     -144.72     50.47                                   
REMARK 500    LYS B  32       78.16     52.07                                   
REMARK 500    VAL B 196      -79.01   -115.74                                   
REMARK 500                                                                      
REMARK 500 REMARK: NULL                                                         
REMARK 620                                                                      
REMARK 620 METAL COORDINATION                                                   
REMARK 620 (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN IDENTIFIER;               
REMARK 620 SSEQ=SEQUENCE NUMBER; I=INSERTION CODE):                             
REMARK 620                                                                      
REMARK 620 COORDINATION ANGLES FOR:  M RES CSSEQI METAL                         
REMARK 620                               K B 302   K                            
REMARK 620 N RES CSSEQI ATOM                                                    
REMARK 620 1 GLN B 223   O                                                      
REMARK 620 2 VAL B 226   O    82.2                                              
REMARK 620 3 Ala B 221   O    86.3  89.8                                        
REMARK 620 4 HOH B 440   O    76.1 158.4  89.6                                  
REMARK 620 5 HOH B 462   O   130.6  60.1  64.6 137.5                            
REMARK 620 6 HOH B 503   O   144.2 129.8 106.5  70.8  84.3                      
REMARK 620 7 HOH B 493   O   114.5  75.3 151.8 113.0  87.2  68.2                
REMARK 620 N                    1     2     3     4     5     6                 
REMARK 620                                                                      
REMARK 620 COORDINATION ANGLES FOR:  M RES CSSEQI METAL                         
REMARK 620                               K A 301   K                            
REMARK 620 N RES CSSEQI ATOM                                                    
REMARK 620 1 VAL A 226   O                                                      
REMARK 620 2 GLN A 223   O    82.4                                              
REMARK 620 3 Ala A 221   O    94.1  82.8                                        
REMARK 620 4 HOH A 519   O    75.0 121.7 150.4                                  
REMARK 620 N                    1     2     3                                   
REMARK 620                                                                      
REMARK 620 COORDINATION ANGLES FOR:  M RES CSSEQI METAL                         
REMARK 620                              NA A 302  NA                            
REMARK 620 N RES CSSEQI ATOM                                                    
REMARK 620 1 LYS A 174   O                                                      
REMARK 620 2 HOH A 497   O   125.8                                              
REMARK 620 N                    1                                               
REMARK 800                                                                      
REMARK 800 SITE                                                                 
REMARK 800 SITE_IDENTIFIER: AC1                                                 
REMARK 800 EVIDENCE_CODE: SOFTWARE                                              
REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE K A 301                   
REMARK 800                                                                      
REMARK 800 SITE_IDENTIFIER: AC2                                                 
REMARK 800 EVIDENCE_CODE: SOFTWARE                                              
REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE NA A 302                  
REMARK 800                                                                      
REMARK 800 SITE_IDENTIFIER: AC3                                                 
REMARK 800 EVIDENCE_CODE: SOFTWARE                                              
REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE BR A 303                  
REMARK 800                                                                      
REMARK 800 SITE_IDENTIFIER: AC4                                                 
REMARK 800 EVIDENCE_CODE: SOFTWARE                                              
REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE PO4 B 301                 
REMARK 800                                                                      
REMARK 800 SITE_IDENTIFIER: AC5                                                 
REMARK 800 EVIDENCE_CODE: SOFTWARE                                              
REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE K B 302                   
REMARK 800                                                                      
REMARK 800 SITE_IDENTIFIER: AC6                                                 
REMARK 800 EVIDENCE_CODE: SOFTWARE                                              
REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE BR B 303                  
REMARK 900                                                                      
REMARK 900 RELATED ENTRIES                                                      
REMARK 900 RELATED ID: 4POD   RELATED DB: PDB                                   
REMARK 900 STRUCTURE OF TRIOSEPHOSPHATE ISOMERASE I170V MUTANT HUMAN            
REMARK 900 ENZYME.                                                              
DBREF  4POC A    0   248  UNP    P60174   TPIS_HUMAN      38    286             
DBREF  4POC B    0   248  UNP    P60174   TPIS_HUMAN      38    286             
SEQADV 4POC Gly A   -5  UNP  P60174              EXPRESSION TAG                 
SEQADV 4POC ASP A   -4  UNP  P60174              EXPRESSION TAG                 
SEQADV 4POC Ile A   -3  UNP  P60174              EXPRESSION TAG                 
SEQADV 4POC Thr A   -2  UNP  P60174              EXPRESSION TAG                 
SEQADV 4POC HIS A   -1  UNP  P60174              EXPRESSION TAG                 
SEQADV 4POC Gly B   -5  UNP  P60174              EXPRESSION TAG                 
SEQADV 4POC ASP B   -4  UNP  P60174              EXPRESSION TAG                 
SEQADV 4POC Ile B   -3  UNP  P60174              EXPRESSION TAG                 
SEQADV 4POC Thr B   -2  UNP  P60174              EXPRESSION TAG                 
SEQADV 4POC HIS B   -1  UNP  P60174              EXPRESSION TAG                 
SEQRES   1 A  254  Gly ASP Ile Thr HIS MET Ala PRO SER ARG LYS Phe PHE          
SEQRES   2 A  254  VAL Gly GLY ASN TRP LYS MET ASN Gly ARG LYS GLN SER          
SEQRES   3 A  254  LEU Gly Glu LEU Ile Gly Thr LEU ASN Ala ALA LYS VAL          
SEQRES   4 A  254  PRO Ala ASP Thr Glu VAL VAL CYS Ala PRO PRO Thr Ala          
SEQRES   5 A  254  TYR Ile ASP Phe Ala ARG GLN LYS LEU ASP PRO LYS Ile          
SEQRES   6 A  254  Ala VAL Ala ALA GLN ASN CYS TYR LYS VAL Thr ASN Gly          
SEQRES   7 A  254  Ala Phe Thr Gly Glu Ile SER PRO Gly MET Ile LYS ASP          
SEQRES   8 A  254  CYS Gly Ala Thr TRP VAL VAL LEU Gly HIS SER Glu ARG          
SEQRES   9 A  254  ARG HIS VAL Phe Gly Glu SER ASP Glu LEU Ile Gly GLN          
SEQRES  10 A  254  LYS VAL Ala HIS Ala LEU Ala Glu Gly LEU Gly VAL Ile          
SEQRES  11 A  254  Ala CYS Ile Gly Glu LYS LEU ASP Glu ARG Glu Ala Gly          
SEQRES  12 A  254  Ile Thr Glu LYS VAL VAL Phe Glu GLN Thr LYS VAL Ile          
SEQRES  13 A  254  Ala ASP ASN VAL LYS ASP TRP SER LYS VAL VAL LEU Ala          
SEQRES  14 A  254  TYR Glu PRO VAL TRP Ala Ile Gly Thr Gly LYS Thr Ala          
SEQRES  15 A  254  Thr PRO GLN GLN Ala GLN Glu VAL HIS Glu LYS LEU ARG          
SEQRES  16 A  254  Gly TRP LEU LYS SER ASN VAL SER ASP Ala VAL Ala GLN          
SEQRES  17 A  254  SER Thr ARG Ile ILE TYR Gly GLY SER VAL Thr Gly Ala          
SEQRES  18 A  254  Thr CYS LYS Glu LEU Ala SER GLN PRO ASP VAL ASP Gly          
SEQRES  19 A  254  Phe LEU VAL Gly GLY Ala SER LEU LYS PRO Glu Phe VAL          
SEQRES  20 A  254  ASP Ile ILE ASN Ala LYS GLN                                  
SEQRES   1 B  254  Gly ASP Ile Thr HIS MET Ala PRO SER ARG LYS Phe PHE          
SEQRES   2 B  254  VAL Gly GLY ASN TRP LYS MET ASN Gly ARG LYS GLN SER          
SEQRES   3 B  254  LEU Gly Glu LEU Ile Gly Thr LEU ASN Ala ALA LYS VAL          
SEQRES   4 B  254  PRO Ala ASP Thr Glu VAL VAL CYS Ala PRO PRO Thr Ala          
SEQRES   5 B  254  TYR Ile ASP Phe Ala ARG GLN LYS LEU ASP PRO LYS Ile          
SEQRES   6 B  254  Ala VAL Ala ALA GLN ASN CYS TYR LYS VAL Thr ASN Gly          
SEQRES   7 B  254  Ala Phe Thr Gly Glu Ile SER PRO Gly MET Ile LYS ASP          
SEQRES   8 B  254  CYS Gly Ala Thr TRP VAL VAL LEU Gly HIS SER Glu ARG          
SEQRES   9 B  254  ARG HIS VAL Phe Gly Glu SER ASP Glu LEU Ile Gly GLN          
SEQRES  10 B  254  LYS VAL Ala HIS Ala LEU Ala Glu Gly LEU Gly VAL Ile          
SEQRES  11 B  254  Ala CYS Ile Gly Glu LYS LEU ASP Glu ARG Glu Ala Gly          
SEQRES  12 B  254  Ile Thr Glu LYS VAL VAL Phe Glu GLN Thr LYS VAL Ile          
SEQRES  13 B  254  Ala ASP ASN VAL LYS ASP TRP SER LYS VAL VAL LEU Ala          
SEQRES  14 B  254  TYR Glu PRO VAL TRP Ala Ile Gly Thr Gly LYS Thr Ala          
SEQRES  15 B  254  Thr PRO GLN GLN Ala GLN Glu VAL HIS Glu LYS LEU ARG          
SEQRES  16 B  254  Gly TRP LEU LYS SER ASN VAL SER ASP Ala VAL Ala GLN          
SEQRES  17 B  254  SER Thr ARG Ile ILE TYR Gly GLY SER VAL Thr Gly Ala          
SEQRES  18 B  254  Thr CYS LYS Glu LEU Ala SER GLN PRO ASP VAL ASP Gly          
SEQRES  19 B  254  Phe LEU VAL Gly GLY Ala SER LEU LYS PRO Glu Phe VAL          
SEQRES  20 B  254  ASP Ile ILE ASN Ala LYS GLN                                  
HET      K  A 301       1                                                       
HET     NA  A 302       1                                                       
HET     BR  A 303       1                                                       
HET    PO4  B 301       5                                                       
HET      K  B 302       1                                                       
HET     BR  B 303       1                                                       
HETNAM       K POTASSIUM ION                                                    
HETNAM      NA SODIUM ION                                                       
HETNAM      BR BROMIDE ION                                                      
HETNAM     PO4 PHOSPHATE ION                                                    
FORMUL   3    K    2(K 1+)                                                      
FORMUL   4   NA    NA 1+                                                        
FORMUL   5   BR    2(BR 1-)                                                     
FORMUL   6  PO4    O4 P 3-                                                      
FORMUL   9  HOH   *327(H2 O)                                                    
HELIX    1   1 ARG A   17  Ala A   31  1                                  15    
HELIX    2   2 PRO A   44  Ala A   46  5                                   3    
HELIX    3   3 TYR A   47  LEU A   55  1                                   9    
HELIX    4   4 SER A   79  CYS A   86  1                                   8    
HELIX    5   5 HIS A   95  VAL A  101  1                                   7    
HELIX    6   6 SER A  105  Glu A  119  1                                  15    
HELIX    7   7 LYS A  130  Ala A  136  1                                   7    
HELIX    8   8 Ile A  138  ASP A  152  1                                  15    
HELIX    9   9 ASP A  156  SER A  158  5                                   3    
HELIX   10  10 PRO A  166  Ile A  170  5                                   5    
HELIX   11  11 Thr A  177  VAL A  196  1                                  20    
HELIX   12  12 SER A  197  Thr A  204  1                                   8    
HELIX   13  13 Thr A  216  SER A  222  1                                   7    
HELIX   14  14 Gly A  232  PRO A  238  5                                   7    
HELIX   15  15 Glu A  239  ASN A  245  1                                   7    
HELIX   16  16 ARG B   17  Ala B   31  1                                  15    
HELIX   17  17 PRO B   44  Ala B   46  5                                   3    
HELIX   18  18 TYR B   47  LEU B   55  1                                   9    
HELIX   19  19 SER B   79  CYS B   86  1                                   8    
HELIX   20  20 HIS B   95  VAL B  101  1                                   7    
HELIX   21  21 SER B  105  Glu B  119  1                                  15    
HELIX   22  22 LYS B  130  Gly B  137  1                                   8    
HELIX   23  23 Ile B  138  ASP B  152  1                                  15    
HELIX   24  24 ASP B  156  SER B  158  5                                   3    
HELIX   25  25 PRO B  166  Ile B  170  5                                   5    
HELIX   26  26 Thr B  177  VAL B  196  1                                  20    
HELIX   27  27 SER B  197  Thr B  204  1                                   8    
HELIX   28  28 Thr B  216  SER B  222  1                                   7    
HELIX   29  29 Gly B  232  LYS B  237  5                                   6    
HELIX   30  30 PRO B  238  ASN B  245  1                                   8    
SHEET    1   A 9 Phe A   6  ASN A  11  0                                        
SHEET    2   A 9 Thr A  37  Ala A  42  1  O  Ala A  42   N  Gly A  10           
SHEET    3   A 9 Ala A  60  Ala A  63  1  O  Ala A  62   N  CYS A  41           
SHEET    4   A 9 TRP A  90  LEU A  93  1  O  VAL A  92   N  Ala A  63           
SHEET    5   A 9 Gly A 122  Ile A 127  1  O  CYS A 126   N  LEU A  93           
SHEET    6   A 9 VAL A 160  TYR A 164  1  O  VAL A 161   N  Ala A 125           
SHEET    7   A 9 Ile A 206  TYR A 208  1  O  Ile A 207   N  LEU A 162           
SHEET    8   A 9 Gly A 228  VAL A 231  1  O  Gly A 228   N  TYR A 208           
SHEET    9   A 9 Phe A   6  ASN A  11  1  N  Gly A   9   O  VAL A 231           
SHEET    1   B 9 Phe B   6  ASN B  11  0                                        
SHEET    2   B 9 Thr B  37  Ala B  42  1  O  Ala B  42   N  Gly B  10           
SHEET    3   B 9 Ala B  60  Ala B  63  1  O  Ala B  62   N  CYS B  41           
SHEET    4   B 9 TRP B  90  LEU B  93  1  O  VAL B  92   N  Ala B  63           
SHEET    5   B 9 Gly B 122  Ile B 127  1  O  CYS B 126   N  LEU B  93           
SHEET    6   B 9 VAL B 160  TYR B 164  1  O  Ala B 163   N  Ala B 125           
SHEET    7   B 9 Ile B 206  Gly B 209  1  O  Ile B 207   N  TYR B 164           
SHEET    8   B 9 Gly B 228  VAL B 231  1  O  Gly B 228   N  TYR B 208           
SHEET    9   B 9 Phe B   6  ASN B  11  1  N  Gly B   9   O  VAL B 231           
LINK         O   GLN B 223                 K     K B 302     1555   1555  2.61  
LINK         O   VAL A 226                 K     K A 301     1555   1555  2.62  
LINK         O   GLN A 223                 K     K A 301     1555   1555  2.64  
LINK         O   VAL B 226                 K     K B 302     1555   1555  2.71  
LINK         O   Ala A 221                 K     K A 301     1555   1555  2.82  
LINK         O   Ala B 221                 K     K B 302     1555   1555  2.83  
LINK         O   LYS A 174                NA    NA A 302     1555   1555  2.93  
LINK         K     K B 302                 O   HOH B 440     1555   1555  2.74  
LINK         K     K B 302                 O   HOH B 462     1555   1555  2.77  
LINK         K     K B 302                 O   HOH B 503     1555   1555  2.79  
LINK         K     K B 302                 O   HOH B 493     1555   1555  2.85  
LINK         K     K A 301                 O   HOH A 519     1555   1555  2.92  
LINK        NA    NA A 302                 O   HOH A 497     1555   1555  3.14  
SITE     1 AC1  4 Ala A 221  GLN A 223  VAL A 226  HOH A 519                    
SITE     1 AC2  6 Ala A 169  Gly A 171  Thr A 172  Gly A 173                    
SITE     2 AC2  6 LYS A 174  HOH B 444                                          
SITE     1 AC3  3 ASN A  71  Thr A 175  HOH A 557                               
SITE     1 AC4 13 LYS B  13  Ala B 169  Ile B 170  Gly B 171                    
SITE     2 AC4 13 Gly B 210  SER B 211  Gly B 232  Gly B 233                    
SITE     3 AC4 13  BR B 303  HOH B 409  HOH B 411  HOH B 477                    
SITE     4 AC4 13 HOH B 481                                                     
SITE     1 AC5  7 Ala B 221  GLN B 223  VAL B 226  HOH B 440                    
SITE     2 AC5  7 HOH B 462  HOH B 493  HOH B 503                               
SITE     1 AC6  4 ASN B  11  HIS B  95  Glu B 165  PO4 B 301                    
CRYST1   47.920   48.851   93.966  90.00 103.66  90.00 P 1 21 1      4          
ORIGX1      1.000000  0.000000  0.000000        0.00000                         
ORIGX2      0.000000  1.000000  0.000000        0.00000                         
ORIGX3      0.000000  0.000000  1.000000        0.00000                         
SCALE1      0.020868  0.000000  0.005072        0.00000                         
SCALE2      0.000000  0.020470  0.000000        0.00000                         
SCALE3      0.000000  0.000000  0.010952        0.00000                         
ATOM   1073  N   Ala A  73       9.852   9.388  28.674  1.00 26.13           N  
ANISOU 1073  N   Ala A  73     3187   3304   3438    390   -712   -860       N  
ATOM   1074  CA  Ala A  73      10.962   9.067  29.557  1.00 25.93           C  
ANISOU 1074  CA  Ala A  73     2965   3468   3419    300   -938   -709       C  
ATOM   1075  C   Ala A  73      11.092   7.560  29.682  1.00 27.19           C  
ANISOU 1075  C   Ala A  73     3217   3692   3420    450   -771   -594       C  
ATOM   1076  O   Ala A  73      10.783   6.978  30.727  1.00 29.60           O  
ANISOU 1076  O   Ala A  73     3747   3949   3550    530   -699   -804       O  
ATOM   1077  CB  Ala A  73      10.741   9.706  30.925  1.00 26.85           C  
ANISOU 1077  CB  Ala A  73     2975   3659   3569    115  -1011   -679       C  
ATOM   1078  H   Ala A  73       9.118   9.542  29.094  1.00 31.36           H  
ATOM   1079  HA  Ala A  73      11.785   9.418  29.183  1.00 31.11           H  
ATOM   1080  HB1 Ala A  73      11.488   9.482  31.501  1.00 32.22           H  
ATOM   1081  HB2 Ala A  73      10.680  10.668  30.818  1.00 32.22           H  
ATOM   1082  HB3 Ala A  73       9.916   9.362  31.303  1.00 32.22           H  
ATOM   1083  N   Phe A  74      11.541   6.937  28.596  1.00 26.02           N  
ANISOU 1083  N   Phe A  74     3013   3549   3325    220   -783   -519       N  
ATOM   1084  CA  Phe A  74      11.684   5.491  28.529  1.00 23.58           C  
ANISOU 1084  CA  Phe A  74     2538   3259   3163    400   -909   -535       C  
ATOM   1085  C   Phe A  74      13.062   5.158  27.985  1.00 22.97           C  
ANISOU 1085  C   Phe A  74     2654   3206   2866    446  -1019   -492       C  
ATOM   1086  O   Phe A  74      13.235   4.752  26.840  1.00 22.87           O  
ANISOU 1086  O   Phe A  74     2461   3347   2880    464  -1067   -622       O  
ATOM   1087  CB  Phe A  74      10.536   4.898  27.709  1.00 23.69           C  
ANISOU 1087  CB  Phe A  74     2385   3266   3349     44   -925   -612       C  
ATOM   1088  CG  Phe A  74       9.192   5.238  28.276  1.00 26.78           C  
ANISOU 1088  CG  Phe A  74     2950   3708   3519   -253   -681   -667       C  
ATOM   1089  CD1 Phe A  74       8.744   4.627  29.435  1.00 28.21           C  
ANISOU 1089  CD1 Phe A  74     2891   4302   3527   -118   -736   -454       C  
ATOM   1090  CD2 Phe A  74       8.401   6.207  27.690  1.00 27.10           C  
ANISOU 1090  CD2 Phe A  74     2757   3955   3586   -200   -833   -682       C  
ATOM   1091  CE1 Phe A  74       7.522   4.958  29.978  1.00 28.51           C  
ANISOU 1091  CE1 Phe A  74     2755   4529   3549   -356   -867   -479       C  
ATOM   1092  CE2 Phe A  74       7.179   6.543  28.235  1.00 26.68           C  
ANISOU 1092  CE2 Phe A  74     2453   4046   3637   -223   -889   -579       C  
ATOM   1093  CZ  Phe A  74       6.745   5.920  29.383  1.00 27.06           C  
ANISOU 1093  CZ  Phe A  74     2345   4262   3675   -367   -908   -636       C  
ATOM   1094  H   Phe A  74      11.773   7.340  27.873  1.00 31.23           H  
ATOM   1095  HA  Phe A  74      11.625   5.128  29.426  1.00 28.30           H  
ATOM   1096  HB2 Phe A  74      10.579   5.247  26.805  1.00 28.42           H  
ATOM   1097  HB3 Phe A  74      10.622   3.932  27.697  1.00 28.42           H  
ATOM   1098  HD1 Phe A  74       9.270   3.980  29.846  1.00 33.86           H  
ATOM   1099  HD2 Phe A  74       8.694   6.635  26.918  1.00 32.53           H  
ATOM   1100  HE1 Phe A  74       7.227   4.536  30.753  1.00 34.22           H  
ATOM   1101  HE2 Phe A  74       6.650   7.190  27.828  1.00 32.01           H  
ATOM   1102  HZ  Phe A  74       5.917   6.138  29.745  1.00 32.47           H  
ATOM   1103  N   Thr A  75      14.045   5.370  28.854  1.00 21.37           N  
ANISOU 1103  N   Thr A  75     2545   2922   2651    219  -1096   -676       N  
ATOM   1104  CA  Thr A  75      15.434   5.077  28.564  1.00 21.59           C  
ANISOU 1104  CA  Thr A  75     2697   3009   2499     40  -1041   -447       C  
ATOM   1105  C   Thr A  75      15.532   3.673  27.995  1.00 21.56           C  
ANISOU 1105  C   Thr A  75     3028   2693   2469    -68   -922   -280       C  
ATOM   1106  O   Thr A  75      14.984   2.730  28.570  1.00 22.34           O  
ANISOU 1106  O   Thr A  75     3232   2894   2362   -134  -1005      6       O  
ATOM   1107  CB  Thr A  75      16.261   5.184  29.859  1.00 22.01           C  
ANISOU 1107  CB  Thr A  75     2774   3073   2515    -75   -940   -516       C  
ATOM   1108  OG1 Thr A  75      16.178   6.525  30.361  1.00 22.88           O  
ANISOU 1108  OG1 Thr A  75     2669   3181   2844    -67   -999   -341       O  
ATOM   1109  CG2 Thr A  75      17.711   4.840  29.627  1.00 23.32           C  
ANISOU 1109  CG2 Thr A  75     3048   3329   2483      3   -620   -527       C  
ATOM   1110  H   Thr A  75      13.923   5.693  29.641  1.00 25.64           H  
ATOM   1111  HA  Thr A  75      15.779   5.708  27.913  1.00 25.91           H  
ATOM   1112  HB  Thr A  75      15.904   4.571  30.520  1.00 26.41           H  
ATOM   1113  HG1 Thr A  75      15.378   6.721  30.525  1.00 27.46           H  
ATOM   1114 HG21 Thr A  75      18.206   4.916  30.458  1.00 27.98           H  
ATOM   1115 HG22 Thr A  75      17.787   3.932  29.297  1.00 27.98           H  
ATOM   1116 HG23 Thr A  75      18.095   5.447  28.975  1.00 27.98           H  
ATOM   1117  N   Gly A  76      16.183   3.555  26.839  1.00 20.82           N  
ANISOU 1117  N   Gly A  76     2890   2528   2493   -283   -948   -386       N  
ATOM   1118  CA  Gly A  76      16.417   2.267  26.209  1.00 21.52           C  
ANISOU 1118  CA  Gly A  76     2887   2716   2573   -135   -945   -319       C  
ATOM   1119  C   Gly A  76      15.364   1.827  25.199  1.00 20.76           C  
ANISOU 1119  C   Gly A  76     2473   2795   2620   -136   -996   -443       C  
ATOM   1120  O   Gly A  76      15.506   0.767  24.588  1.00 21.88           O  
ANISOU 1120  O   Gly A  76     2407   3199   2709    319   -996   -497       O  
ATOM   1121  H   Gly A  76      16.502   4.220  26.397  1.00 24.98           H  
ATOM   1122  HA2 Gly A  76      17.272   2.294  25.752  1.00 25.82           H  
ATOM   1123  HA3 Gly A  76      16.468   1.587  26.898  1.00 25.82           H  
ATOM   1124  N   Glu A  77      14.314   2.627  25.017  1.00 20.38           N  
ANISOU 1124  N   Glu A  77     2274   2789   2680    -29   -932   -370       N  
ATOM   1125  CA  Glu A  77      13.250   2.296  24.058  1.00 20.85           C  
ANISOU 1125  CA  Glu A  77     2267   2914   2741    308  -1152   -278       C  
ATOM   1126  C   Glu A  77      13.303   3.173  22.816  1.00 21.36           C  
ANISOU 1126  C   Glu A  77     2608   2655   2851    356  -1201   -248       C  
ATOM   1127  O   Glu A  77      13.968   4.209  22.802  1.00 21.76           O  
ANISOU 1127  O   Glu A  77     2862   2658   2750     27  -1328   -180       O  
ATOM   1128  CB  Glu A  77      11.868   2.470  24.701  1.00 23.44           C  
ANISOU 1128  CB  Glu A  77     2710   3344   2854    261   -929   -377       C  
ATOM   1129  CG  Glu A  77      11.656   1.620  25.940  1.00 24.69           C  
ANISOU 1129  CG  Glu A  77     2767   3655   2958    130   -836   -321       C  
ATOM   1130  CD  Glu A  77      11.643   0.145  25.620  1.00 25.80           C  
ANISOU 1130  CD  Glu A  77     3031   3781   2990   -166   -930   -288       C  
ATOM   1131  OE1 Glu A  77      10.846  -0.259  24.739  1.00 27.63           O  
ANISOU 1131  OE1 Glu A  77     3536   3791   3173   -241   -928   -217       O  
ATOM   1132  OE2 Glu A  77      12.427  -0.609  26.239  1.00 25.79           O  
ANISOU 1132  OE2 Glu A  77     3043   3826   2932      0   -929   -323       O  
ATOM   1133  H   Glu A  77      14.191   3.369  25.434  1.00 24.45           H  
ATOM   1134  HA  Glu A  77      13.344   1.370  23.782  1.00 25.02           H  
ATOM   1135  HB2 Glu A  77      11.758   3.399  24.956  1.00 28.13           H  
ATOM   1136  HB3 Glu A  77      11.190   2.224  24.053  1.00 28.13           H  
ATOM   1137  HG2 Glu A  77      12.376   1.786  26.567  1.00 29.62           H  
ATOM   1138  HG3 Glu A  77      10.803   1.849  26.341  1.00 29.62           H  
ATOM   1139  N   Ile A  78      12.606   2.735  21.770  1.00 20.06           N  
ANISOU 1139  N   Ile A  78     2488   2313   2820    599  -1487   -220       N  
ATOM   1140  CA  Ile A  78      12.384   3.553  20.579  1.00 20.59           C  
ANISOU 1140  CA  Ile A  78     2436   2494   2895    297  -1526   -263       C  
ATOM   1141  C   Ile A  78      10.902   3.545  20.209  1.00 22.20           C  
ANISOU 1141  C   Ile A  78     2601   2819   3016    402  -1357   -200       C  
ATOM   1142  O   Ile A  78      10.180   2.608  20.539  1.00 23.69           O  
ANISOU 1142  O   Ile A  78     2405   3421   3177    660  -1010   -242       O  
ATOM   1143  CB  Ile A  78      13.219   3.068  19.369  1.00 23.39           C  
ANISOU 1143  CB  Ile A  78     2838   2890   3159    600   -954   -141       C  
ATOM   1144  CG1 Ile A  78      12.812   1.654  18.927  1.00 23.81           C  
ANISOU 1144  CG1 Ile A  78     2646   3145   3256    564  -1052   -550       C  
ATOM   1145  CG2 Ile A  78      14.721   3.134  19.690  1.00 22.80           C  
ANISOU 1145  CG2 Ile A  78     2890   2822   2950    835  -1169    224       C  
ATOM   1146  CD1 Ile A  78      13.346   1.276  17.558  1.00 24.35           C  
ANISOU 1146  CD1 Ile A  78     2771   3262   3219    453  -1424   -569       C  
ATOM   1147  H   Ile A  78      12.247   1.955  21.724  1.00 24.07           H  
ATOM   1148  HA  Ile A  78      12.640   4.468  20.773  1.00 24.71           H  
ATOM   1149  HB  Ile A  78      13.047   3.671  18.629  1.00 28.07           H  
ATOM   1150 HG12 Ile A  78      13.156   1.013  19.569  1.00 28.58           H  
ATOM   1151 HG13 Ile A  78      11.844   1.601  18.895  1.00 28.58           H  
ATOM   1152 HG21 Ile A  78      15.222   2.826  18.919  1.00 27.36           H  
ATOM   1153 HG22 Ile A  78      14.959   4.052  19.895  1.00 27.36           H  
ATOM   1154 HG23 Ile A  78      14.905   2.565  20.454  1.00 27.36           H  
ATOM   1155 HD11 Ile A  78      13.053   0.377  17.343  1.00 29.22           H  
ATOM   1156 HD12 Ile A  78      13.002   1.903  16.902  1.00 29.22           H  
ATOM   1157 HD13 Ile A  78      14.315   1.314  17.577  1.00 29.22           H                                                                                                                     
MASTER      310    0    6   30   18    0   11    6 4056    2   22   40          
END  */