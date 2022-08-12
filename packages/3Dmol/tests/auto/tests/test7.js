
/*

@div
<div style="width: 400px; height: 400px; position: relative;"
               class='viewer_3Dmoljs'
               data-element='moldata_pdb'
               data-backgroundcolor='0xffffff'
               data-select1='chain:A' data-style1='cartoon:style=trace' data-surface1='opacity:.7;color:red'
               data-select2='chain:B' data-style2='cartoon:thickness=0.1'
               data-type="pdb"></div>

*/

/* 

  @data   <textarea style="display: none;" id="moldata_pdb">
			HEADER    UNKNOWN FUNCTION                        17-SEP-08   3EIT              
TITLE     THE 2.6 ANGSTROM CRYSTAL STRUCTURE OF CHBP, THE CIF HOMOLOGUE FROM    
TITLE    2 BURKHOLDERIA PSEUDOMALLEI                                            
COMPND    MOL_ID: 1;                                                            
COMPND   2 MOLECULE: PUTATIVE ATP/GTP BINDING PROTEIN;                          
COMPND   3 CHAIN: A, B;                                                         
COMPND   4 FRAGMENT: RESIDUES 48-328;                                           
COMPND   5 SYNONYM: THE EPEC EFFECTOR CIF HOMOLOGUE;                            
COMPND   6 ENGINEERED: YES                                                      
SOURCE    MOL_ID: 1;                                                            
SOURCE   2 ORGANISM_SCIENTIFIC: BURKHOLDERIA PSEUDOMALLEI;                      
SOURCE   3 ORGANISM_TAXID: 28450;                                               
SOURCE   4 STRAIN: K96243;                                                      
SOURCE   5 GENE: YP_111397;                                                     
SOURCE   6 EXPRESSION_SYSTEM: ESCHERICHIA COLI;                                 
SOURCE   7 EXPRESSION_SYSTEM_TAXID: 562;                                        
SOURCE   8 EXPRESSION_SYSTEM_STRAIN: BL21(DE3);                                 
SOURCE   9 EXPRESSION_SYSTEM_VECTOR_TYPE: PLASMID;                              
SOURCE  10 EXPRESSION_SYSTEM_PLASMID: PGEX-6P-2                                 
KEYWDS    PAPAIN-LIKE FOLD,APAIN SUPERFAMILY, HYDROLYTIC ENZYME, UNKNOWN        
KEYWDS   2 FUNCTION                                                             
EXPDTA    X-RAY DIFFRACTION                                                     
AUTHOR    Q.YAO,Y.ZHU,F.SHAO                                                    
REVDAT   3   20-JUL-11 3EIT    1       REMARK                                   
REVDAT   2   07-APR-09 3EIT    1       JRNL                                     
REVDAT   1   03-FEB-09 3EIT    0                                                
JRNL        AUTH   Q.YAO,J.CUI,Y.ZHU,G.WANG,L.HU,C.LONG,R.CAO,X.LIU,N.HUANG,    
JRNL        AUTH 2 S.CHEN,L.LIU,F.SHAO                                          
JRNL        TITL   A BACTERIAL TYPE III EFFECTOR FAMILY USES THE PAPAIN-LIKE    
JRNL        TITL 2 HYDROLYTIC ACTIVITY TO ARREST THE HOST CELL CYCLE            
JRNL        REF    PROC.NATL.ACAD.SCI.USA        V. 106  3716 2009              
JRNL        REFN                   ISSN 0027-8424                               
JRNL        PMID   19225106                                                     
JRNL        DOI    10.1073/PNAS.0900212106                                      
REMARK   2                                                                      
REMARK   2 RESOLUTION.    2.56 ANGSTROMS.                                       
REMARK   3                                                                      
REMARK   3 REFINEMENT.                                                          
REMARK   3   PROGRAM     : REFMAC 5.2.0019                                      
REMARK   3   AUTHORS     : MURSHUDOV,VAGIN,DODSON                               
REMARK   3                                                                      
REMARK   3    REFINEMENT TARGET : MAXIMUM LIKELIHOOD                            
REMARK   3                                                                      
REMARK   3  DATA USED IN REFINEMENT.                                            
REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : 2.56                           
REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : 49.15                          
REMARK   3   DATA CUTOFF            (SIGMA(F)) : 0.000                          
REMARK   3   COMPLETENESS FOR RANGE        (%) : 98.1                           
REMARK   3   NUMBER OF REFLECTIONS             : 16115                          
REMARK   3                                                                      
REMARK   3  FIT TO DATA USED IN REFINEMENT.                                     
REMARK   3   CROSS-VALIDATION METHOD          : THROUGHOUT                      
REMARK   3   FREE R VALUE TEST SET SELECTION  : RANDOM                          
REMARK   3   R VALUE     (WORKING + TEST SET) : 0.210                           
REMARK   3   R VALUE            (WORKING SET) : 0.203                           
REMARK   3   FREE R VALUE                     : 0.267                           
REMARK   3   FREE R VALUE TEST SET SIZE   (%) : 10.100                          
REMARK   3   FREE R VALUE TEST SET COUNT      : 1635                            
REMARK   3                                                                      
REMARK   3  FIT IN THE HIGHEST RESOLUTION BIN.                                  
REMARK   3   TOTAL NUMBER OF BINS USED           : 20                           
REMARK   3   BIN RESOLUTION RANGE HIGH       (A) : 2.56                         
REMARK   3   BIN RESOLUTION RANGE LOW        (A) : 2.62                         
REMARK   3   REFLECTION IN BIN     (WORKING SET) : 804                          
REMARK   3   BIN COMPLETENESS (WORKING+TEST) (%) : 75.57                        
REMARK   3   BIN R VALUE           (WORKING SET) : 0.2650                       
REMARK   3   BIN FREE R VALUE SET COUNT          : 87                           
REMARK   3   BIN FREE R VALUE                    : 0.3900                       
REMARK   3                                                                      
REMARK   3  NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT.                    
REMARK   3   PROTEIN ATOMS            : 3853                                    
REMARK   3   NUCLEIC ACID ATOMS       : 0                                       
REMARK   3   HETEROGEN ATOMS          : 0                                       
REMARK   3   SOLVENT ATOMS            : 70                                      
REMARK   3                                                                      
REMARK   3  B VALUES.                                                           
REMARK   3   B VALUE TYPE : UNVERIFIED                                          
REMARK   3   FROM WILSON PLOT           (A**2) : 51.50                          
REMARK   3   MEAN B VALUE      (OVERALL, A**2) : 30.50                          
REMARK   3   OVERALL ANISOTROPIC B VALUE.                                       
REMARK   3    B11 (A**2) : 4.02000                                              
REMARK   3    B22 (A**2) : -1.77000                                             
REMARK   3    B33 (A**2) : -2.25000                                             
REMARK   3    B12 (A**2) : 0.00000                                              
REMARK   3    B13 (A**2) : 0.00000                                              
REMARK   3    B23 (A**2) : 0.00000                                              
REMARK   3                                                                      
REMARK   3  ESTIMATED OVERALL COORDINATE ERROR.                                 
REMARK   3   ESU BASED ON R VALUE                            (A): NULL          
REMARK   3   ESU BASED ON FREE R VALUE                       (A): 0.371         
REMARK   3   ESU BASED ON MAXIMUM LIKELIHOOD                 (A): 0.247         
REMARK   3   ESU FOR B VALUES BASED ON MAXIMUM LIKELIHOOD (A**2): 23.894        
REMARK   3                                                                      
REMARK   3 CORRELATION COEFFICIENTS.                                            
REMARK   3   CORRELATION COEFFICIENT FO-FC      : 0.940                         
REMARK   3   CORRELATION COEFFICIENT FO-FC FREE : 0.903                         
REMARK   3                                                                      
REMARK   3  RMS DEVIATIONS FROM IDEAL VALUES        COUNT    RMS    WEIGHT      
REMARK   3   BOND LENGTHS REFINED ATOMS        (A):  3922 ; 0.008 ; 0.022       
REMARK   3   BOND LENGTHS OTHERS               (A):  NULL ;  NULL ;  NULL       
REMARK   3   BOND ANGLES REFINED ATOMS   (DEGREES):  5310 ; 1.110 ; 1.969       
REMARK   3   BOND ANGLES OTHERS          (DEGREES):  NULL ;  NULL ;  NULL       
REMARK   3   TORSION ANGLES, PERIOD 1    (DEGREES):   484 ; 5.580 ; 5.000       
REMARK   3   TORSION ANGLES, PERIOD 2    (DEGREES):   186 ;35.811 ;24.355       
REMARK   3   TORSION ANGLES, PERIOD 3    (DEGREES):   694 ;17.946 ;15.000       
REMARK   3   TORSION ANGLES, PERIOD 4    (DEGREES):    28 ;15.970 ;15.000       
REMARK   3   CHIRAL-CENTER RESTRAINTS       (A**3):   600 ; 0.071 ; 0.200       
REMARK   3   GENERAL PLANES REFINED ATOMS      (A):  2948 ; 0.003 ; 0.020       
REMARK   3   GENERAL PLANES OTHERS             (A):  NULL ;  NULL ;  NULL       
REMARK   3   NON-BONDED CONTACTS REFINED ATOMS (A):  1699 ; 0.201 ; 0.200       
REMARK   3   NON-BONDED CONTACTS OTHERS        (A):  NULL ;  NULL ;  NULL       
REMARK   3   NON-BONDED TORSION REFINED ATOMS  (A):  2720 ; 0.297 ; 0.200       
REMARK   3   NON-BONDED TORSION OTHERS         (A):  NULL ;  NULL ;  NULL       
REMARK   3   H-BOND (X...Y) REFINED ATOMS      (A):   165 ; 0.140 ; 0.200       
REMARK   3   H-BOND (X...Y) OTHERS             (A):  NULL ;  NULL ;  NULL       
REMARK   3   POTENTIAL METAL-ION REFINED ATOMS (A):  NULL ;  NULL ;  NULL       
REMARK   3   POTENTIAL METAL-ION OTHERS        (A):  NULL ;  NULL ;  NULL       
REMARK   3   SYMMETRY VDW REFINED ATOMS        (A):    63 ; 0.212 ; 0.200       
REMARK   3   SYMMETRY VDW OTHERS               (A):  NULL ;  NULL ;  NULL       
REMARK   3   SYMMETRY H-BOND REFINED ATOMS     (A):    12 ; 0.240 ; 0.200       
REMARK   3   SYMMETRY H-BOND OTHERS            (A):  NULL ;  NULL ;  NULL       
REMARK   3   SYMMETRY METAL-ION REFINED ATOMS  (A):  NULL ;  NULL ;  NULL       
REMARK   3   SYMMETRY METAL-ION OTHERS         (A):  NULL ;  NULL ;  NULL       
REMARK   3                                                                      
REMARK   3  ISOTROPIC THERMAL FACTOR RESTRAINTS.     COUNT   RMS    WEIGHT      
REMARK   3   MAIN-CHAIN BOND REFINED ATOMS  (A**2):  2485 ; 0.401 ; 1.500       
REMARK   3   MAIN-CHAIN BOND OTHER ATOMS    (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   MAIN-CHAIN ANGLE REFINED ATOMS (A**2):  3916 ; 0.652 ; 2.000       
REMARK   3   SIDE-CHAIN BOND REFINED ATOMS  (A**2):  1587 ; 1.283 ; 3.000       
REMARK   3   SIDE-CHAIN ANGLE REFINED ATOMS (A**2):  1394 ; 1.919 ; 4.500       
REMARK   3                                                                      
REMARK   3 ANISOTROPIC THERMAL FACTOR RESTRAINTS.    COUNT   RMS   WEIGHT       
REMARK   3   RIGID-BOND RESTRAINTS          (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   SPHERICITY; FREE ATOMS         (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   SPHERICITY; BONDED ATOMS       (A**2):  NULL ;  NULL ;  NULL       
REMARK   3                                                                      
REMARK   3  NCS RESTRAINTS STATISTICS                                           
REMARK   3   NUMBER OF DIFFERENT NCS GROUPS : NULL                              
REMARK   3                                                                      
REMARK   3  TLS DETAILS                                                         
REMARK   3   NUMBER OF TLS GROUPS  : 14                                         
REMARK   3                                                                      
REMARK   3   TLS GROUP : 1                                                      
REMARK   3    NUMBER OF COMPONENTS GROUP : 1                                    
REMARK   3    COMPONENTS        C SSSEQI   TO  C SSSEQI                         
REMARK   3    RESIDUE RANGE :   A    79        A   119                          
REMARK   3    ORIGIN FOR THE GROUP (A):  39.0661  34.5399  64.0379              
REMARK   3    T TENSOR                                                          
REMARK   3      T11:   0.0794 T22:   0.0671                                     
REMARK   3      T33:   0.0212 T12:  -0.0056                                     
REMARK   3      T13:   0.0132 T23:   0.0088                                     
REMARK   3    L TENSOR                                                          
REMARK   3      L11:   3.2977 L22:   5.5308                                     
REMARK   3      L33:   2.7537 L12:  -1.5386                                     
REMARK   3      L13:   0.2530 L23:  -1.4130                                     
REMARK   3    S TENSOR                                                          
REMARK   3      S11:  -0.0558 S12:  -0.0984 S13:   0.2448                       
REMARK   3      S21:   0.1365 S22:   0.1019 S23:   0.0745                       
REMARK   3      S31:  -0.0560 S32:  -0.0586 S33:  -0.0461                       
REMARK   3                                                                      
REMARK   3   TLS GROUP : 2                                                      
REMARK   3    NUMBER OF COMPONENTS GROUP : 1                                    
REMARK   3    COMPONENTS        C SSSEQI   TO  C SSSEQI                         
REMARK   3    RESIDUE RANGE :   A   120        A   151                          
REMARK   3    ORIGIN FOR THE GROUP (A):  21.6302  39.5226  49.9442              
REMARK   3    T TENSOR                                                          
REMARK   3      T11:   0.0308 T22:   0.0697                                     
REMARK   3      T33:   0.0946 T12:   0.0421                                     
REMARK   3      T13:   0.0745 T23:   0.0546                                     
REMARK   3    L TENSOR                                                          
REMARK   3      L11:   4.6614 L22:   8.6481                                     
REMARK   3      L33:   6.8617 L12:  -1.4222                                     
REMARK   3      L13:   0.4337 L23:  -1.9140                                     
REMARK   3    S TENSOR                                                          
REMARK   3      S11:   0.2995 S12:   0.1447 S13:  -0.1036                       
REMARK   3      S21:  -0.2408 S22:  -0.0529 S23:   0.2177                       
REMARK   3      S31:   0.2160 S32:  -0.1487 S33:  -0.2465                       
REMARK   3                                                                      
REMARK   3   TLS GROUP : 3                                                      
REMARK   3    NUMBER OF COMPONENTS GROUP : 1                                    
REMARK   3    COMPONENTS        C SSSEQI   TO  C SSSEQI                         
REMARK   3    RESIDUE RANGE :   A   152        A   181                          
REMARK   3    ORIGIN FOR THE GROUP (A):  44.7432  40.1613  56.7824              
REMARK   3    T TENSOR                                                          
REMARK   3      T11:   0.1442 T22:  -0.0017                                     
REMARK   3      T33:   0.0369 T12:   0.0685                                     
REMARK   3      T13:   0.0444 T23:   0.0406                                     
REMARK   3    L TENSOR                                                          
REMARK   3      L11:   3.0196 L22:   2.7470                                     
REMARK   3      L33:   4.2011 L12:  -0.9439                                     
REMARK   3      L13:  -0.1787 L23:   0.1720                                     
REMARK   3    S TENSOR                                                          
REMARK   3      S11:   0.0143 S12:   0.0637 S13:   0.1779                       
REMARK   3      S21:  -0.0075 S22:   0.0657 S23:  -0.2546                       
REMARK   3      S31:  -0.3579 S32:  -0.0090 S33:  -0.0800                       
REMARK   3                                                                      
REMARK   3   TLS GROUP : 4                                                      
REMARK   3    NUMBER OF COMPONENTS GROUP : 1                                    
REMARK   3    COMPONENTS        C SSSEQI   TO  C SSSEQI                         
REMARK   3    RESIDUE RANGE :   A   182        A   197                          
REMARK   3    ORIGIN FOR THE GROUP (A):  63.8946  37.1862  55.7769              
REMARK   3    T TENSOR                                                          
REMARK   3      T11:  -0.0193 T22:   0.0538                                     
REMARK   3      T33:   0.1473 T12:  -0.0173                                     
REMARK   3      T13:   0.0198 T23:   0.0419                                     
REMARK   3    L TENSOR                                                          
REMARK   3      L11:   8.4475 L22:   8.3084                                     
REMARK   3      L33:   3.4474 L12:  -3.5701                                     
REMARK   3      L13:  -4.1923 L23:   2.0605                                     
REMARK   3    S TENSOR                                                          
REMARK   3      S11:   0.2208 S12:   0.1949 S13:   0.4872                       
REMARK   3      S21:   0.2496 S22:  -0.3423 S23:  -0.3717                       
REMARK   3      S31:   0.2979 S32:   0.1324 S33:   0.1215                       
REMARK   3                                                                      
REMARK   3   TLS GROUP : 5                                                      
REMARK   3    NUMBER OF COMPONENTS GROUP : 1                                    
REMARK   3    COMPONENTS        C SSSEQI   TO  C SSSEQI                         
REMARK   3    RESIDUE RANGE :   A   198        A   249                          
REMARK   3    ORIGIN FOR THE GROUP (A):  49.4501  28.7886  54.4754              
REMARK   3    T TENSOR                                                          
REMARK   3      T11:   0.1122 T22:   0.0534                                     
REMARK   3      T33:   0.0505 T12:   0.0587                                     
REMARK   3      T13:   0.0206 T23:   0.0094                                     
REMARK   3    L TENSOR                                                          
REMARK   3      L11:   2.4896 L22:   2.2456                                     
REMARK   3      L33:   2.3118 L12:  -0.1493                                     
REMARK   3      L13:   0.8982 L23:  -0.6426                                     
REMARK   3    S TENSOR                                                          
REMARK   3      S11:   0.0422 S12:   0.0750 S13:   0.0891                       
REMARK   3      S21:  -0.1357 S22:  -0.1709 S23:  -0.0562                       
REMARK   3      S31:   0.1361 S32:   0.1858 S33:   0.1287                       
REMARK   3                                                                      
REMARK   3   TLS GROUP : 6                                                      
REMARK   3    NUMBER OF COMPONENTS GROUP : 1                                    
REMARK   3    COMPONENTS        C SSSEQI   TO  C SSSEQI                         
REMARK   3    RESIDUE RANGE :   A   250        A   272                          
REMARK   3    ORIGIN FOR THE GROUP (A):  58.6240  24.4573  50.6790              
REMARK   3    T TENSOR                                                          
REMARK   3      T11:   0.2287 T22:   0.0160                                     
REMARK   3      T33:   0.1458 T12:   0.0606                                     
REMARK   3      T13:   0.1139 T23:  -0.0396                                     
REMARK   3    L TENSOR                                                          
REMARK   3      L11:   4.4178 L22:   4.2101                                     
REMARK   3      L33:   3.4691 L12:  -3.3311                                     
REMARK   3      L13:   0.1149 L23:   2.3396                                     
REMARK   3    S TENSOR                                                          
REMARK   3      S11:   0.1490 S12:   0.6437 S13:  -0.0160                       
REMARK   3      S21:  -0.3999 S22:  -0.2526 S23:  -0.8757                       
REMARK   3      S31:  -0.2981 S32:   0.0200 S33:   0.1036                       
REMARK   3                                                                      
REMARK   3   TLS GROUP : 7                                                      
REMARK   3    NUMBER OF COMPONENTS GROUP : 1                                    
REMARK   3    COMPONENTS        C SSSEQI   TO  C SSSEQI                         
REMARK   3    RESIDUE RANGE :   A   273        A   318                          
REMARK   3    ORIGIN FOR THE GROUP (A):  52.6986  26.0919  50.0159              
REMARK   3    T TENSOR                                                          
REMARK   3      T11:   0.1251 T22:   0.0334                                     
REMARK   3      T33:  -0.0263 T12:   0.0034                                     
REMARK   3      T13:  -0.0150 T23:   0.0415                                     
REMARK   3    L TENSOR                                                          
REMARK   3      L11:   2.6356 L22:   3.9414                                     
REMARK   3      L33:   1.9622 L12:  -0.9446                                     
REMARK   3      L13:   0.0225 L23:   0.4123                                     
REMARK   3    S TENSOR                                                          
REMARK   3      S11:   0.0766 S12:   0.5409 S13:  -0.1555                       
REMARK   3      S21:  -0.6620 S22:  -0.1379 S23:  -0.0223                       
REMARK   3      S31:   0.0364 S32:   0.1636 S33:   0.0613                       
REMARK   3                                                                      
REMARK   3   TLS GROUP : 8                                                      
REMARK   3    NUMBER OF COMPONENTS GROUP : 1                                    
REMARK   3    COMPONENTS        C SSSEQI   TO  C SSSEQI                         
REMARK   3    RESIDUE RANGE :   A   319        A   327                          
REMARK   3    ORIGIN FOR THE GROUP (A):  37.6125  36.9913  70.2708              
REMARK   3    T TENSOR                                                          
REMARK   3      T11:   0.4260 T22:   0.0872                                     
REMARK   3      T33:   0.0257 T12:   0.0687                                     
REMARK   3      T13:   0.0088 T23:  -0.0805                                     
REMARK   3    L TENSOR                                                          
REMARK   3      L11:  19.9181 L22:   9.7321                                     
REMARK   3      L33:  19.2786 L12:  -1.6263                                     
REMARK   3      L13:  -8.9989 L23: -11.3497                                     
REMARK   3    S TENSOR                                                          
REMARK   3      S11:   0.4350 S12:  -0.9308 S13:   0.6723                       
REMARK   3      S21:   2.1256 S22:  -0.3085 S23:   0.1661                       
REMARK   3      S31:  -0.7779 S32:  -0.2601 S33:  -0.1265                       
REMARK   3                                                                      
REMARK   3   TLS GROUP : 9                                                      
REMARK   3    NUMBER OF COMPONENTS GROUP : 1                                    
REMARK   3    COMPONENTS        C SSSEQI   TO  C SSSEQI                         
REMARK   3    RESIDUE RANGE :   B    80        B   119                          
REMARK   3    ORIGIN FOR THE GROUP (A):  36.8273  42.4579  29.1226              
REMARK   3    T TENSOR                                                          
REMARK   3      T11:   0.2364 T22:   0.0445                                     
REMARK   3      T33:  -0.0485 T12:   0.0246                                     
REMARK   3      T13:   0.0077 T23:  -0.0137                                     
REMARK   3    L TENSOR                                                          
REMARK   3      L11:   2.9211 L22:   5.4818                                     
REMARK   3      L33:   3.7365 L12:  -0.3321                                     
REMARK   3      L13:   1.0160 L23:  -2.7960                                     
REMARK   3    S TENSOR                                                          
REMARK   3      S11:  -0.1014 S12:  -0.1248 S13:  -0.0284                       
REMARK   3      S21:  -0.0771 S22:   0.1275 S23:   0.0500                       
REMARK   3      S31:   0.1330 S32:  -0.1000 S33:  -0.0261                       
REMARK   3                                                                      
REMARK   3   TLS GROUP : 10                                                     
REMARK   3    NUMBER OF COMPONENTS GROUP : 1                                    
REMARK   3    COMPONENTS        C SSSEQI   TO  C SSSEQI                         
REMARK   3    RESIDUE RANGE :   B   120        B   151                          
REMARK   3    ORIGIN FOR THE GROUP (A):  32.9336  24.2227  40.3415              
REMARK   3    T TENSOR                                                          
REMARK   3      T11:   0.1489 T22:   0.0045                                     
REMARK   3      T33:   0.0357 T12:   0.0283                                     
REMARK   3      T13:  -0.0566 T23:  -0.0162                                     
REMARK   3    L TENSOR                                                          
REMARK   3      L11:   2.9986 L22:   6.3269                                     
REMARK   3      L33:   5.8231 L12:   0.1783                                     
REMARK   3      L13:  -0.6550 L23:  -2.5160                                     
REMARK   3    S TENSOR                                                          
REMARK   3      S11:  -0.0547 S12:   0.1434 S13:  -0.0474                       
REMARK   3      S21:  -0.3374 S22:   0.0109 S23:   0.1810                       
REMARK   3      S31:  -0.0271 S32:  -0.2365 S33:   0.0438                       
REMARK   3                                                                      
REMARK   3   TLS GROUP : 11                                                     
REMARK   3    NUMBER OF COMPONENTS GROUP : 1                                    
REMARK   3    COMPONENTS        C SSSEQI   TO  C SSSEQI                         
REMARK   3    RESIDUE RANGE :   B   152        B   185                          
REMARK   3    ORIGIN FOR THE GROUP (A):  44.3012  47.7087  38.5836              
REMARK   3    T TENSOR                                                          
REMARK   3      T11:   0.0848 T22:   0.0381                                     
REMARK   3      T33:  -0.0308 T12:   0.0239                                     
REMARK   3      T13:  -0.0135 T23:   0.0544                                     
REMARK   3    L TENSOR                                                          
REMARK   3      L11:   3.0609 L22:   4.7012                                     
REMARK   3      L33:   6.9136 L12:   0.1110                                     
REMARK   3      L13:  -1.0651 L23:   0.5354                                     
REMARK   3    S TENSOR                                                          
REMARK   3      S11:   0.0683 S12:  -0.2971 S13:   0.0552                       
REMARK   3      S21:   0.3146 S22:  -0.0465 S23:  -0.4995                       
REMARK   3      S31:  -0.4018 S32:   0.8816 S33:  -0.0218                       
REMARK   3                                                                      
REMARK   3   TLS GROUP : 12                                                     
REMARK   3    NUMBER OF COMPONENTS GROUP : 1                                    
REMARK   3    COMPONENTS        C SSSEQI   TO  C SSSEQI                         
REMARK   3    RESIDUE RANGE :   B   186        B   291                          
REMARK   3    ORIGIN FOR THE GROUP (A):  34.1202  59.6853  41.2713              
REMARK   3    T TENSOR                                                          
REMARK   3      T11:   0.2376 T22:   0.0225                                     
REMARK   3      T33:   0.0936 T12:  -0.0010                                     
REMARK   3      T13:   0.0133 T23:   0.0284                                     
REMARK   3    L TENSOR                                                          
REMARK   3      L11:   2.4176 L22:   1.2992                                     
REMARK   3      L33:   3.9444 L12:  -1.1018                                     
REMARK   3      L13:  -1.2942 L23:   1.9457                                     
REMARK   3    S TENSOR                                                          
REMARK   3      S11:  -0.0186 S12:  -0.0609 S13:   0.4230                       
REMARK   3      S21:  -0.0435 S22:   0.1180 S23:  -0.0752                       
REMARK   3      S31:  -0.4999 S32:   0.3411 S33:  -0.0994                       
REMARK   3                                                                      
REMARK   3   TLS GROUP : 13                                                     
REMARK   3    NUMBER OF COMPONENTS GROUP : 1                                    
REMARK   3    COMPONENTS        C SSSEQI   TO  C SSSEQI                         
REMARK   3    RESIDUE RANGE :   B   292        B   306                          
REMARK   3    ORIGIN FOR THE GROUP (A):  30.7993  56.4838  53.9843              
REMARK   3    T TENSOR                                                          
REMARK   3      T11:   0.3786 T22:   0.0170                                     
REMARK   3      T33:   0.0842 T12:   0.0411                                     
REMARK   3      T13:   0.0270 T23:  -0.0218                                     
REMARK   3    L TENSOR                                                          
REMARK   3      L11:  14.9130 L22:   4.1489                                     
REMARK   3      L33:  10.0751 L12:   0.2380                                     
REMARK   3      L13:  -0.7658 L23:  -1.3452                                     
REMARK   3    S TENSOR                                                          
REMARK   3      S11:  -0.8503 S12:  -1.2381 S13:  -0.4136                       
REMARK   3      S21:   0.4657 S22:   0.4840 S23:  -0.0189                       
REMARK   3      S31:  -0.2612 S32:  -0.0787 S33:   0.3663                       
REMARK   3                                                                      
REMARK   3   TLS GROUP : 14                                                     
REMARK   3    NUMBER OF COMPONENTS GROUP : 1                                    
REMARK   3    COMPONENTS        C SSSEQI   TO  C SSSEQI                         
REMARK   3    RESIDUE RANGE :   B   307        B   327                          
REMARK   3    ORIGIN FOR THE GROUP (A):  43.1045  49.1468  26.6277              
REMARK   3    T TENSOR                                                          
REMARK   3      T11:   0.1391 T22:   0.0352                                     
REMARK   3      T33:  -0.0202 T12:   0.0615                                     
REMARK   3      T13:   0.0401 T23:   0.1109                                     
REMARK   3    L TENSOR                                                          
REMARK   3      L11:   5.1865 L22:   4.7139                                     
REMARK   3      L33:   5.9898 L12:   2.4428                                     
REMARK   3      L13:  -0.4874 L23:   0.5424                                     
REMARK   3    S TENSOR                                                          
REMARK   3      S11:  -0.0718 S12:   0.5755 S13:   0.1069                       
REMARK   3      S21:  -0.4480 S22:   0.2929 S23:  -0.3048                       
REMARK   3      S31:  -0.0793 S32:   0.6485 S33:  -0.2211                       
REMARK   3                                                                      
REMARK   3  BULK SOLVENT MODELLING.                                             
REMARK   3   METHOD USED : MASK                                                 
REMARK   3   PARAMETERS FOR MASK CALCULATION                                    
REMARK   3   VDW PROBE RADIUS   : 1.20                                          
REMARK   3   ION PROBE RADIUS   : 0.80                                          
REMARK   3   SHRINKAGE RADIUS   : 0.80                                          
REMARK   3                                                                      
REMARK   3  OTHER REFINEMENT REMARKS: HYDROGENS HAVE BEEN ADDED IN THE RIDING   
REMARK   3  POSITIONS                                                           
REMARK   4                                                                      
REMARK   4 3EIT COMPLIES WITH FORMAT V. 3.30, 13-JUL-11                         
REMARK 100                                                                      
REMARK 100 THIS ENTRY HAS BEEN PROCESSED BY PDBJ ON 29-SEP-08.                  
REMARK 100 THE RCSB ID CODE IS RCSB049375.                                      
REMARK 200                                                                      
REMARK 200 EXPERIMENTAL DETAILS                                                 
REMARK 200  EXPERIMENT TYPE                : X-RAY DIFFRACTION                  
REMARK 200  DATE OF DATA COLLECTION        : 03-MAY-08                          
REMARK 200  TEMPERATURE           (KELVIN) : 100                                
REMARK 200  PH                             : 6.0                                
REMARK 200  NUMBER OF CRYSTALS USED        : 1                                  
REMARK 200                                                                      
REMARK 200  SYNCHROTRON              (Y/N) : Y                                  
REMARK 200  RADIATION SOURCE               : PHOTON FACTORY                     
REMARK 200  BEAMLINE                       : BL-5A                              
REMARK 200  X-RAY GENERATOR MODEL          : NULL                               
REMARK 200  MONOCHROMATIC OR LAUE    (M/L) : M                                  
REMARK 200  WAVELENGTH OR RANGE        (A) : 0.97912                            
REMARK 200  MONOCHROMATOR                  : NULL                               
REMARK 200  OPTICS                         : NULL                               
REMARK 200                                                                      
REMARK 200  DETECTOR TYPE                  : CCD                                
REMARK 200  DETECTOR MANUFACTURER          : RIGAKU                             
REMARK 200  INTENSITY-INTEGRATION SOFTWARE : HKL-2000                           
REMARK 200  DATA SCALING SOFTWARE          : HKL-2000                           
REMARK 200                                                                      
REMARK 200  NUMBER OF UNIQUE REFLECTIONS   : 16433                              
REMARK 200  RESOLUTION RANGE HIGH      (A) : 2.560                              
REMARK 200  RESOLUTION RANGE LOW       (A) : 50.000                             
REMARK 200  REJECTION CRITERIA  (SIGMA(I)) : 1.000                              
REMARK 200                                                                      
REMARK 200 OVERALL.                                                             
REMARK 200  COMPLETENESS FOR RANGE     (%) : 100.0                              
REMARK 200  DATA REDUNDANCY                : 7.400                              
REMARK 200  R MERGE                    (I) : 0.07400                            
REMARK 200  R SYM                      (I) : 0.07400                            
REMARK 200  <I/SIGMA(I)> FOR THE DATA SET  : 26.9000                            
REMARK 200                                                                      
REMARK 200 IN THE HIGHEST RESOLUTION SHELL.                                     
REMARK 200  HIGHEST RESOLUTION SHELL, RANGE HIGH (A) : 2.56                     
REMARK 200  HIGHEST RESOLUTION SHELL, RANGE LOW  (A) : 2.62                     
REMARK 200  COMPLETENESS FOR SHELL     (%) : 75.0                               
REMARK 200  DATA REDUNDANCY IN SHELL       : 6.70                               
REMARK 200  R MERGE FOR SHELL          (I) : 0.27800                            
REMARK 200  R SYM FOR SHELL            (I) : 0.27800                            
REMARK 200  <I/SIGMA(I)> FOR SHELL         : 6.400                              
REMARK 200                                                                      
REMARK 200 DIFFRACTION PROTOCOL: SINGLE WAVELENGTH                              
REMARK 200 METHOD USED TO DETERMINE THE STRUCTURE: SAD                          
REMARK 200 SOFTWARE USED: SOLVE                                                 
REMARK 200 STARTING MODEL: NULL                                                 
REMARK 200                                                                      
REMARK 200 REMARK: NULL                                                         
REMARK 280                                                                      
REMARK 280 CRYSTAL                                                              
REMARK 280 SOLVENT CONTENT, VS   (%): 36.28                                     
REMARK 280 MATTHEWS COEFFICIENT, VM (ANGSTROMS**3/DA): 1.93                     
REMARK 280                                                                      
REMARK 280 CRYSTALLIZATION CONDITIONS: 32% PEG1000, 100MM SODIUM CACODYLATE,    
REMARK 280  5% GLYCEROL, PH 6.0, VAPOR DIFFUSION, HANGING DROP, TEMPERATURE     
REMARK 280  293K                                                                
REMARK 290                                                                      
REMARK 290 CRYSTALLOGRAPHIC SYMMETRY                                            
REMARK 290 SYMMETRY OPERATORS FOR SPACE GROUP: P 21 21 21                       
REMARK 290                                                                      
REMARK 290      SYMOP   SYMMETRY                                                
REMARK 290     NNNMMM   OPERATOR                                                
REMARK 290       1555   X,Y,Z                                                   
REMARK 290       2555   -X+1/2,-Y,Z+1/2                                         
REMARK 290       3555   -X,Y+1/2,-Z+1/2                                         
REMARK 290       4555   X+1/2,-Y+1/2,-Z                                         
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
REMARK 290   SMTRY1   2 -1.000000  0.000000  0.000000       27.18100            
REMARK 290   SMTRY2   2  0.000000 -1.000000  0.000000        0.00000            
REMARK 290   SMTRY3   2  0.000000  0.000000  1.000000       57.51850            
REMARK 290   SMTRY1   3 -1.000000  0.000000  0.000000        0.00000            
REMARK 290   SMTRY2   3  0.000000  1.000000  0.000000       39.00150            
REMARK 290   SMTRY3   3  0.000000  0.000000 -1.000000       57.51850            
REMARK 290   SMTRY1   4  1.000000  0.000000  0.000000       27.18100            
REMARK 290   SMTRY2   4  0.000000 -1.000000  0.000000       39.00150            
REMARK 290   SMTRY3   4  0.000000  0.000000 -1.000000        0.00000            
REMARK 290                                                                      
REMARK 290 REMARK: NULL                                                         
REMARK 300                                                                      
REMARK 300 BIOMOLECULE: 1, 2                                                    
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
REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: MONOMERIC                         
REMARK 350 APPLY THE FOLLOWING TO CHAINS: A                                     
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000            
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000            
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000            
REMARK 350                                                                      
REMARK 350 BIOMOLECULE: 2                                                       
REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: MONOMERIC                         
REMARK 350 APPLY THE FOLLOWING TO CHAINS: B                                     
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
REMARK 465     GLY A    48                                                      
REMARK 465     LEU A    49                                                      
REMARK 465     PRO A    50                                                      
REMARK 465     ALA A    51                                                      
REMARK 465     ARG A    52                                                      
REMARK 465     SER A    53                                                      
REMARK 465     SER A    54                                                      
REMARK 465     SER A    55                                                      
REMARK 465     ILE A    56                                                      
REMARK 465     SER A    57                                                      
REMARK 465     ASN A    58                                                      
REMARK 465     THR A    59                                                      
REMARK 465     ASN A    60                                                      
REMARK 465     ARG A    61                                                      
REMARK 465     THR A    62                                                      
REMARK 465     GLY A    63                                                      
REMARK 465     GLU A    64                                                      
REMARK 465     ASN A    65                                                      
REMARK 465     PRO A    66                                                      
REMARK 465     MSE A    67                                                      
REMARK 465     ILE A    68                                                      
REMARK 465     THR A    69                                                      
REMARK 465     PRO A    70                                                      
REMARK 465     ILE A    71                                                      
REMARK 465     ILE A    72                                                      
REMARK 465     SER A    73                                                      
REMARK 465     SER A    74                                                      
REMARK 465     ASN A    75                                                      
REMARK 465     LEU A    76                                                      
REMARK 465     GLY A    77                                                      
REMARK 465     LEU A    78                                                      
REMARK 465     LEU A   220                                                      
REMARK 465     THR A   221                                                      
REMARK 465     GLN A   222                                                      
REMARK 465     GLY A   223                                                      
REMARK 465     PRO A   224                                                      
REMARK 465     GLY A   328                                                      
REMARK 465     GLY B    48                                                      
REMARK 465     LEU B    49                                                      
REMARK 465     PRO B    50                                                      
REMARK 465     ALA B    51                                                      
REMARK 465     ARG B    52                                                      
REMARK 465     SER B    53                                                      
REMARK 465     SER B    54                                                      
REMARK 465     SER B    55                                                      
REMARK 465     ILE B    56                                                      
REMARK 465     SER B    57                                                      
REMARK 465     ASN B    58                                                      
REMARK 465     THR B    59                                                      
REMARK 465     ASN B    60                                                      
REMARK 465     ARG B    61                                                      
REMARK 465     THR B    62                                                      
REMARK 465     GLY B    63                                                      
REMARK 465     GLU B    64                                                      
REMARK 465     ASN B    65                                                      
REMARK 465     PRO B    66                                                      
REMARK 465     MSE B    67                                                      
REMARK 465     ILE B    68                                                      
REMARK 465     THR B    69                                                      
REMARK 465     PRO B    70                                                      
REMARK 465     ILE B    71                                                      
REMARK 465     ILE B    72                                                      
REMARK 465     SER B    73                                                      
REMARK 465     SER B    74                                                      
REMARK 465     ASN B    75                                                      
REMARK 465     LEU B    76                                                      
REMARK 465     GLY B    77                                                      
REMARK 465     LEU B    78                                                      
REMARK 465     LYS B    79                                                      
REMARK 465     LEU B   220                                                      
REMARK 465     THR B   221                                                      
REMARK 465     GLN B   222                                                      
REMARK 465     GLY B   223                                                      
REMARK 465     GLY B   328                                                      
REMARK 470                                                                      
REMARK 470 MISSING ATOM                                                         
REMARK 470 THE FOLLOWING RESIDUES HAVE MISSING ATOMS (M=MODEL NUMBER;           
REMARK 470 RES=RESIDUE NAME; C=CHAIN IDENTIFIER; SSEQ=SEQUENCE NUMBER;          
REMARK 470 I=INSERTION CODE):                                                   
REMARK 470   M RES CSSEQI  ATOMS                                                
REMARK 470     HIS B  80    CG   ND1  CD2  CE1  NE2                             
REMARK 470     ARG B  81    CG   CD   NE   CZ   NH1  NH2                        
REMARK 470     ARG B 109    CG   CD   NE   CZ   NH1  NH2                        
REMARK 470     LYS B 267    CG   CD   CE   NZ                                   
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
REMARK 500    LYS A 170      -78.46   -114.92                                   
REMARK 500    LYS B 170      -86.13   -110.64                                   
REMARK 500    THR B 177       22.24    -78.79                                   
REMARK 500    SER B 178       -4.26   -147.27                                   
REMARK 500    GLU B 183      -26.90   -145.60                                   
REMARK 500                                                                      
REMARK 500 REMARK: NULL                                                         
REMARK 900                                                                      
REMARK 900 RELATED ENTRIES                                                      
REMARK 900 RELATED ID: 3EIR   RELATED DB: PDB                                   
DBREF  3EIT A   48   328  UNP    Q63KH5   Q63KH5_BURPS    48    328             
DBREF  3EIT B   48   328  UNP    Q63KH5   Q63KH5_BURPS    48    328             
SEQRES   1 A  281  GLY LEU PRO ALA ARG SER SER SER ILE SER ASN THR ASN          
SEQRES   2 A  281  ARG THR GLY GLU ASN PRO MSE ILE THR PRO ILE ILE SER          
SEQRES   3 A  281  SER ASN LEU GLY LEU LYS HIS ARG VAL THR LEU ARG LYS          
SEQRES   4 A  281  ALA THR LEU ALA SER LEU MSE GLN SER LEU SER GLY GLU          
SEQRES   5 A  281  SER SER ASN ARG VAL MSE TRP ASN ASP ARG TYR ASP THR          
SEQRES   6 A  281  LEU LEU ILE ALA ARG ASP PRO ARG GLU ILE LYS ASN ALA          
SEQRES   7 A  281  ILE GLU LYS SER VAL THR ASP PHE GLY GLY LEU GLU ASN          
SEQRES   8 A  281  TYR LYS GLU LEU THR GLY GLY ALA ASP PRO PHE ALA LEU          
SEQRES   9 A  281  MSE THR PRO VAL OCS GLY LEU SER ALA ASN ASN ILE PHE          
SEQRES  10 A  281  LYS LEU MSE THR GLU LYS ASP VAL PRO ILE ASP PRO THR          
SEQRES  11 A  281  SER ILE GLU TYR LEU GLU ASN THR SER PHE ALA GLU HIS          
SEQRES  12 A  281  VAL ASN THR LEU ASP SER HIS LYS ASN TYR VAL VAL ILE          
SEQRES  13 A  281  VAL ASN ASP GLY ARG LEU GLY HIS LYS PHE LEU ILE ASP          
SEQRES  14 A  281  LEU PRO ALA LEU THR GLN GLY PRO ARG THR ALA TYR ILE          
SEQRES  15 A  281  ILE GLN SER ASP LEU GLY GLY GLY ALA LEU PRO ALA VAL          
SEQRES  16 A  281  ARG VAL GLU ASP TRP ILE SER ARG ARG GLY SER ASP PRO          
SEQRES  17 A  281  VAL SER LEU ASP GLU LEU ASN GLN LEU LEU SER LYS ASP          
SEQRES  18 A  281  PHE SER LYS MSE PRO ASP ASP VAL GLN THR ARG LEU LEU          
SEQRES  19 A  281  ALA SER ILE LEU GLN ILE ASP LYS ASP PRO HIS LYS VAL          
SEQRES  20 A  281  ASP ILE LYS LYS LEU HIS LEU ASP GLY LYS LEU ARG PHE          
SEQRES  21 A  281  ALA SER HIS GLU TYR ASP PHE ARG GLN PHE GLN ARG ASN          
SEQRES  22 A  281  ALA GLN TYR VAL ALA GLY LEU GLY                              
SEQRES   1 B  281  GLY LEU PRO ALA ARG SER SER SER ILE SER ASN THR ASN          
SEQRES   2 B  281  ARG THR GLY GLU ASN PRO MSE ILE THR PRO ILE ILE SER          
SEQRES   3 B  281  SER ASN LEU GLY LEU LYS HIS ARG VAL THR LEU ARG LYS          
SEQRES   4 B  281  ALA THR LEU ALA SER LEU MSE GLN SER LEU SER GLY GLU          
SEQRES   5 B  281  SER SER ASN ARG VAL MSE TRP ASN ASP ARG TYR ASP THR          
SEQRES   6 B  281  LEU LEU ILE ALA ARG ASP PRO ARG GLU ILE LYS ASN ALA          
SEQRES   7 B  281  ILE GLU LYS SER VAL THR ASP PHE GLY GLY LEU GLU ASN          
SEQRES   8 B  281  TYR LYS GLU LEU THR GLY GLY ALA ASP PRO PHE ALA LEU          
SEQRES   9 B  281  MSE THR PRO VAL OCS GLY LEU SER ALA ASN ASN ILE PHE          
SEQRES  10 B  281  LYS LEU MSE THR GLU LYS ASP VAL PRO ILE ASP PRO THR          
SEQRES  11 B  281  SER ILE GLU TYR LEU GLU ASN THR SER PHE ALA GLU HIS          
SEQRES  12 B  281  VAL ASN THR LEU ASP SER HIS LYS ASN TYR VAL VAL ILE          
SEQRES  13 B  281  VAL ASN ASP GLY ARG LEU GLY HIS LYS PHE LEU ILE ASP          
SEQRES  14 B  281  LEU PRO ALA LEU THR GLN GLY PRO ARG THR ALA TYR ILE          
SEQRES  15 B  281  ILE GLN SER ASP LEU GLY GLY GLY ALA LEU PRO ALA VAL          
SEQRES  16 B  281  ARG VAL GLU ASP TRP ILE SER ARG ARG GLY SER ASP PRO          
SEQRES  17 B  281  VAL SER LEU ASP GLU LEU ASN GLN LEU LEU SER LYS ASP          
SEQRES  18 B  281  PHE SER LYS MSE PRO ASP ASP VAL GLN THR ARG LEU LEU          
SEQRES  19 B  281  ALA SER ILE LEU GLN ILE ASP LYS ASP PRO HIS LYS VAL          
SEQRES  20 B  281  ASP ILE LYS LYS LEU HIS LEU ASP GLY LYS LEU ARG PHE          
SEQRES  21 B  281  ALA SER HIS GLU TYR ASP PHE ARG GLN PHE GLN ARG ASN          
SEQRES  22 B  281  ALA GLN TYR VAL ALA GLY LEU GLY                              
MODRES 3EIT MSE A   93  MET  SELENOMETHIONINE                                   
MODRES 3EIT MSE A  105  MET  SELENOMETHIONINE                                   
MODRES 3EIT MSE A  152  MET  SELENOMETHIONINE                                   
MODRES 3EIT OCS A  156  CYS  CYSTEINESULFONIC ACID                              
MODRES 3EIT MSE A  167  MET  SELENOMETHIONINE                                   
MODRES 3EIT MSE A  272  MET  SELENOMETHIONINE                                   
MODRES 3EIT MSE B   93  MET  SELENOMETHIONINE                                   
MODRES 3EIT MSE B  105  MET  SELENOMETHIONINE                                   
MODRES 3EIT MSE B  152  MET  SELENOMETHIONINE                                   
MODRES 3EIT OCS B  156  CYS  CYSTEINESULFONIC ACID                              
MODRES 3EIT MSE B  167  MET  SELENOMETHIONINE                                   
MODRES 3EIT MSE B  272  MET  SELENOMETHIONINE                                   
HET    MSE  A  93       8                                                       
HET    MSE  A 105       8                                                       
HET    MSE  A 152       8                                                       
HET    OCS  A 156       9                                                       
HET    MSE  A 167       8                                                       
HET    MSE  A 272       8                                                       
HET    MSE  B  93       8                                                       
HET    MSE  B 105       8                                                       
HET    MSE  B 152       8                                                       
HET    OCS  B 156       9                                                       
HET    MSE  B 167       8                                                       
HET    MSE  B 272       8                                                       
HETNAM     MSE SELENOMETHIONINE                                                 
HETNAM     OCS CYSTEINESULFONIC ACID                                            
FORMUL   1  MSE    10(C5 H11 N O2 SE)                                           
FORMUL   1  OCS    2(C3 H7 N O5 S)                                              
FORMUL   3  HOH   *70(H2 O)                                                     
HELIX    1   1 HIS A   80  LEU A   96  1                                  17    
HELIX    2   2 GLY A   98  ASN A  107  1                                  10    
HELIX    3   3 LEU A  113  ARG A  117  5                                   5    
HELIX    4   4 ASP A  118  PHE A  133  1                                  16    
HELIX    5   5 GLY A  135  GLY A  144  1                                  10    
HELIX    6   6 VAL A  155  GLU A  169  1                                  15    
HELIX    7   7 SER A  186  LEU A  194  1                                   9    
HELIX    8   8 ARG A  243  GLY A  252  1                                  10    
HELIX    9   9 SER A  257  LEU A  265  1                                   9    
HELIX   10  10 SER A  266  MSE A  272  5                                   7    
HELIX   11  11 PRO A  273  GLN A  286  1                                  14    
HELIX   12  12 ASP A  290  VAL A  294  5                                   5    
HELIX   13  13 ASP A  295  LEU A  299  5                                   5    
HELIX   14  14 ASP A  313  GLY A  326  1                                  14    
HELIX   15  15 HIS B   80  LEU B   96  1                                  17    
HELIX   16  16 GLY B   98  TRP B  106  1                                   9    
HELIX   17  17 LEU B  113  ARG B  117  5                                   5    
HELIX   18  18 ASP B  118  PHE B  133  1                                  16    
HELIX   19  19 GLY B  135  THR B  143  1                                   9    
HELIX   20  20 VAL B  155  GLU B  169  1                                  15    
HELIX   21  21 ASP B  175  ILE B  179  5                                   5    
HELIX   22  22 SER B  186  ASN B  192  1                                   7    
HELIX   23  23 ARG B  243  GLY B  252  1                                  10    
HELIX   24  24 SER B  257  LEU B  265  1                                   9    
HELIX   25  25 ASP B  268  MSE B  272  5                                   5    
HELIX   26  26 PRO B  273  GLN B  286  1                                  14    
HELIX   27  27 ASP B  290  VAL B  294  5                                   5    
HELIX   28  28 ASP B  295  LEU B  299  5                                   5    
HELIX   29  29 ASP B  313  GLY B  326  1                                  14    
SHEET    1   A 4 TYR A 228  ILE A 230  0                                        
SHEET    2   A 4 HIS A 211  LEU A 217 -1  N  LEU A 214   O  ILE A 230           
SHEET    3   A 4 ASN A 199  ASP A 206 -1  N  VAL A 202   O  ILE A 215           
SHEET    4   A 4 LEU A 305  TYR A 312 -1  O  ALA A 308   N  ILE A 203           
SHEET    1   B 4 TYR B 228  ILE B 230  0                                        
SHEET    2   B 4 HIS B 211  LEU B 217 -1  N  ASP B 216   O  TYR B 228           
SHEET    3   B 4 ASN B 199  ASP B 206 -1  N  TYR B 200   O  LEU B 217           
SHEET    4   B 4 LEU B 305  TYR B 312 -1  O  ARG B 306   N  ASN B 205           
LINK         C   LEU A  92                 N   MSE A  93     1555   1555  1.31  
LINK         C   MSE A  93                 N   GLN A  94     1555   1555  1.31  
LINK         C   VAL A 104                 N   MSE A 105     1555   1555  1.33  
LINK         C   MSE A 105                 N   TRP A 106     1555   1555  1.34  
LINK         C   LEU A 151                 N   MSE A 152     1555   1555  1.33  
LINK         C   MSE A 152                 N   THR A 153     1555   1555  1.33  
LINK         C   VAL A 155                 N   OCS A 156     1555   1555  1.33  
LINK         C   OCS A 156                 N   GLY A 157     1555   1555  1.33  
LINK         C   LEU A 166                 N   MSE A 167     1555   1555  1.33  
LINK         C   MSE A 167                 N   THR A 168     1555   1555  1.33  
LINK         C   LYS A 271                 N   MSE A 272     1555   1555  1.33  
LINK         C   MSE A 272                 N   PRO A 273     1555   1555  1.35  
LINK         C   LEU B  92                 N   MSE B  93     1555   1555  1.33  
LINK         C   MSE B  93                 N   GLN B  94     1555   1555  1.34  
LINK         C   VAL B 104                 N   MSE B 105     1555   1555  1.32  
LINK         C   MSE B 105                 N   TRP B 106     1555   1555  1.31  
LINK         C   LEU B 151                 N   MSE B 152     1555   1555  1.33  
LINK         C   MSE B 152                 N   THR B 153     1555   1555  1.33  
LINK         C   VAL B 155                 N   OCS B 156     1555   1555  1.33  
LINK         C   OCS B 156                 N   GLY B 157     1555   1555  1.33  
LINK         C   LEU B 166                 N   MSE B 167     1555   1555  1.33  
LINK         C   MSE B 167                 N   THR B 168     1555   1555  1.32  
LINK         C   LYS B 271                 N   MSE B 272     1555   1555  1.33  
LINK         C   MSE B 272                 N   PRO B 273     1555   1555  1.35  
CRYST1   54.362   78.003  115.037  90.00  90.00  90.00 P 21 21 21    8          
ORIGX1      1.000000  0.000000  0.000000        0.00000                         
ORIGX2      0.000000  1.000000  0.000000        0.00000                         
ORIGX3      0.000000  0.000000  1.000000        0.00000                         
SCALE1      0.018395  0.000000  0.000000        0.00000                         
SCALE2      0.000000  0.012820  0.000000        0.00000                         
SCALE3      0.000000  0.000000  0.008693        0.00000                         
ATOM      1  N   LYS A  79      58.865  26.297  76.719  1.00 31.70           N  
ATOM      2  CA  LYS A  79      57.636  27.111  76.478  1.00 31.83           C  
ATOM      3  C   LYS A  79      56.478  26.289  75.914  1.00 31.75           C  
ATOM      4  O   LYS A  79      56.666  25.178  75.397  1.00 31.76           O  
ATOM      5  CB  LYS A  79      57.936  28.319  75.576  1.00 31.92           C  
ATOM      6  CG  LYS A  79      58.642  29.449  76.309  1.00 32.06           C  
ATOM      7  CD  LYS A  79      58.740  30.702  75.461  1.00 32.18           C  
ATOM      8  CE  LYS A  79      58.611  31.946  76.329  1.00 31.97           C  
ATOM      9  NZ  LYS A  79      57.189  32.174  76.724  1.00 31.68           N  
ATOM     10  N   HIS A  80      55.284  26.862  76.011  1.00 31.41           N  
ATOM     11  CA  HIS A  80      54.059  26.161  75.693  1.00 31.06           C  
ATOM     12  C   HIS A  80      53.283  26.898  74.629  1.00 30.77           C  
ATOM     13  O   HIS A  80      52.064  26.742  74.503  1.00 30.89           O  
ATOM     14  CB  HIS A  80      53.218  26.019  76.956  1.00 31.34           C  
ATOM     15  CG  HIS A  80      53.839  25.137  77.995  1.00 31.69           C  
ATOM     16  ND1 HIS A  80      54.463  25.634  79.120  1.00 31.32           N  
ATOM     17  CD2 HIS A  80      53.948  23.789  78.067  1.00 31.28           C  
ATOM     18  CE1 HIS A  80      54.921  24.631  79.845  1.00 31.22           C  
ATOM     19  NE2 HIS A  80      54.621  23.501  79.229  1.00 31.81           N  
ATOM     20  N   ARG A  81      54.000  27.703  73.856  1.00 30.28           N  
ATOM     21  CA  ARG A  81      53.405  28.430  72.757  1.00 29.83           C  
ATOM     22  C   ARG A  81      52.644  27.504  71.804  1.00 29.11           C  
ATOM     23  O   ARG A  81      51.523  27.806  71.407  1.00 29.22           O  
ATOM     24  CB  ARG A  81      54.481  29.193  71.991  1.00 30.25           C  
ATOM     25  CG  ARG A  81      53.938  30.330  71.119  1.00 32.05           C  
ATOM     26  CD  ARG A  81      54.103  31.672  71.804  1.00 35.30           C  
ATOM     27  NE  ARG A  81      55.514  32.034  71.898  1.00 38.32           N  
ATOM     28  CZ  ARG A  81      55.970  33.208  72.320  1.00 39.58           C  
ATOM     29  NH1 ARG A  81      55.127  34.165  72.702  1.00 39.12           N  
ATOM     30  NH2 ARG A  81      57.280  33.422  72.358  1.00 40.53           N  
ATOM     31  N   VAL A  82      53.245  26.376  71.444  1.00 28.23           N  
ATOM     32  CA  VAL A  82      52.652  25.499  70.437  1.00 27.44           C  
ATOM     33  C   VAL A  82      51.315  24.918  70.905  1.00 26.78           C  
ATOM     34  O   VAL A  82      50.350  24.892  70.133  1.00 26.60           O  
ATOM     35  CB  VAL A  82      53.634  24.394  69.980  1.00 27.54           C  
ATOM     36  CG1 VAL A  82      53.033  23.572  68.848  1.00 27.34           C  
ATOM     37  CG2 VAL A  82      54.963  25.018  69.533  1.00 27.39           C  
ATOM     38  N   THR A  83      51.267  24.475  72.163  1.00 25.89           N  
ATOM     39  CA  THR A  83      50.051  23.933  72.775  1.00 25.30           C  
ATOM     40  C   THR A  83      48.945  24.981  72.787  1.00 25.13           C  
ATOM     41  O   THR A  83      47.787  24.693  72.465  1.00 25.68           O  
ATOM     42  CB  THR A  83      50.301  23.461  74.233  1.00 25.62           C  
ATOM     43  OG1 THR A  83      51.465  22.619  74.288  1.00 26.57           O  
ATOM     44  CG2 THR A  83      49.089  22.694  74.789  1.00 24.54           C  
ATOM     45  N   LEU A  84      49.310  26.204  73.146  1.00 24.39           N  
ATOM     46  CA  LEU A  84      48.354  27.289  73.251  1.00 23.52           C  
ATOM     47  C   LEU A  84      47.867  27.713  71.877  1.00 23.14           C  
ATOM     48  O   LEU A  84      46.687  28.012  71.673  1.00 22.88           O  
ATOM     49  CB  LEU A  84      48.997  28.461  73.994  1.00 23.48           C  
ATOM     50  CG  LEU A  84      49.118  28.235  75.498  1.00 22.62           C  
ATOM     51  CD1 LEU A  84      50.099  29.207  76.134  1.00 21.83           C  
ATOM     52  CD2 LEU A  84      47.743  28.355  76.124  1.00 22.69           C  
ATOM     53  N   ARG A  85      48.800  27.734  70.940  1.00 22.78           N  
ATOM     54  CA  ARG A  85      48.527  28.088  69.558  1.00 22.52           C  
ATOM     55  C   ARG A  85      47.531  27.090  68.949  1.00 22.54           C  
ATOM     56  O   ARG A  85      46.582  27.481  68.262  1.00 22.11           O  
ATOM     57  CB  ARG A  85      49.853  28.091  68.806  1.00 22.43           C  
ATOM     58  CG  ARG A  85      49.749  28.170  67.336  1.00 21.94           C  
ATOM     59  CD  ARG A  85      51.121  28.082  66.686  1.00 19.92           C  
ATOM     60  NE  ARG A  85      50.995  28.504  65.292  1.00 19.97           N  
ATOM     61  CZ  ARG A  85      50.490  27.754  64.313  1.00 18.11           C  
ATOM     62  NH1 ARG A  85      50.088  26.521  64.554  1.00 16.17           N  
ATOM     63  NH2 ARG A  85      50.406  28.240  63.080  1.00 19.08           N  
ATOM     64  N   LYS A  86      47.746  25.806  69.240  1.00 22.79           N  
ATOM     65  CA  LYS A  86      46.852  24.738  68.809  1.00 22.99           C  
ATOM     66  C   LYS A  86      45.452  24.912  69.417  1.00 23.06           C  
ATOM     67  O   LYS A  86      44.454  24.827  68.709  1.00 22.81           O  
ATOM     68  CB  LYS A  86      47.453  23.388  69.186  1.00 23.19           C  
ATOM     69  CG  LYS A  86      46.667  22.163  68.730  1.00 23.83           C  
ATOM     70  CD  LYS A  86      47.383  20.916  69.217  1.00 26.62           C  
ATOM     71  CE  LYS A  86      46.765  19.655  68.650  1.00 30.50           C  
ATOM     72  NZ  LYS A  86      46.815  19.655  67.148  1.00 33.61           N  
ATOM     73  N   ALA A  87      45.394  25.169  70.725  1.00 23.05           N  
ATOM     74  CA  ALA A  87      44.132  25.448  71.406  1.00 23.30           C  
ATOM     75  C   ALA A  87      43.446  26.679  70.829  1.00 23.37           C  
ATOM     76  O   ALA A  87      42.252  26.657  70.542  1.00 23.65           O  
ATOM     77  CB  ALA A  87      44.358  25.614  72.913  1.00 23.32           C  
ATOM     78  N   THR A  88      44.207  27.752  70.658  1.00 23.75           N  
ATOM     79  CA  THR A  88      43.670  28.981  70.090  1.00 24.27           C  
ATOM     80  C   THR A  88      43.109  28.754  68.679  1.00 24.31           C  
ATOM     81  O   THR A  88      41.979  29.150  68.401  1.00 24.74           O  
ATOM     82  CB  THR A  88      44.715  30.119  70.100  1.00 24.34           C  
ATOM     83  OG1 THR A  88      45.163  30.342  71.442  1.00 24.78           O  
ATOM     84  CG2 THR A  88      44.122  31.405  69.565  1.00 24.64           C  
ATOM     85  N   LEU A  89      43.878  28.110  67.805  1.00 24.01           N  
ATOM     86  CA  LEU A  89      43.397  27.836  66.447  1.00 24.12           C  
ATOM     87  C   LEU A  89      42.182  26.899  66.457  1.00 24.30           C  
ATOM     88  O   LEU A  89      41.235  27.097  65.701  1.00 24.08           O  
ATOM     89  CB  LEU A  89      44.519  27.276  65.556  1.00 23.62           C  
ATOM     90  CG  LEU A  89      45.663  28.230  65.199  1.00 23.37           C  
ATOM     91  CD1 LEU A  89      46.864  27.478  64.616  1.00 22.91           C  
ATOM     92  CD2 LEU A  89      45.215  29.352  64.276  1.00 21.52           C  
ATOM     93  N   ALA A  90      42.216  25.889  67.322  1.00 24.87           N  
ATOM     94  CA  ALA A  90      41.112  24.937  67.460  1.00 25.79           C  
ATOM     95  C   ALA A  90      39.831  25.659  67.828  1.00 26.54           C  
ATOM     96  O   ALA A  90      38.859  25.628  67.074  1.00 26.77           O  
ATOM     97  CB  ALA A  90      41.431  23.869  68.519  1.00 25.40           C  
ATOM     98  N   SER A  91      39.846  26.321  68.984  1.00 27.42           N  
ATOM     99  CA  SER A  91      38.660  26.980  69.500  1.00 28.42           C  
ATOM    100  C   SER A  91      38.140  27.983  68.469  1.00 28.41           C  
ATOM    101  O   SER A  91      36.928  28.109  68.282  1.00 29.26           O  
ATOM    102  CB  SER A  91      38.936  27.638  70.858  1.00 28.28           C  
ATOM    103  OG  SER A  91      39.022  29.047  70.752  1.00 30.70           O  
ATOM    104  N   LEU A  92      39.052  28.662  67.781  1.00 28.12           N  
ATOM    105  CA  LEU A  92      38.658  29.574  66.708  1.00 28.41           C  
ATOM    106  C   LEU A  92      37.935  28.837  65.566  1.00 28.10           C  
ATOM    107  O   LEU A  92      36.842  29.244  65.167  1.00 27.93           O  
ATOM    108  CB  LEU A  92      39.855  30.391  66.182  1.00 27.79           C  
ATOM    109  CG  LEU A  92      39.544  31.364  65.025  1.00 29.17           C  
ATOM    110  CD1 LEU A  92      38.669  32.582  65.418  1.00 28.28           C  
ATOM    111  CD2 LEU A  92      40.829  31.818  64.324  1.00 29.02           C  
HETATM  112  N   MSE A  93      38.500  27.756  65.077  1.00 20.00           N  
HETATM  113  CA  MSE A  93      37.902  27.032  63.990  1.00 20.00           C  
HETATM  114  C   MSE A  93      36.575  26.453  64.301  1.00 20.00           C  
HETATM  115  O   MSE A  93      35.741  26.397  63.470  1.00 28.07           O  
HETATM  116  CB  MSE A  93      38.803  25.920  63.483  1.00 20.00           C  
HETATM  117  CG  MSE A  93      40.036  26.350  62.755  1.00 20.00           C  
HETATM  118 SE   MSE A  93      39.890  27.368  61.208  1.00 20.00          SE  
HETATM  119  CE  MSE A  93      39.972  29.069  61.872  1.00 20.00           C  
ATOM    120  N   GLN A  94      36.420  25.958  65.499  1.00 28.49           N  
ATOM    121  CA  GLN A  94      35.166  25.350  65.922  1.00 29.26           C  
ATOM    122  C   GLN A  94      34.013  26.343  65.850  1.00 29.19           C  
ATOM    123  O   GLN A  94      32.886  25.964  65.534  1.00 29.44           O  
ATOM    124  CB  GLN A  94      35.273  24.814  67.352  1.00 29.36           C  
ATOM    125  CG  GLN A  94      36.180  23.601  67.554  1.00 29.68           C  
ATOM    126  CD  GLN A  94      36.497  23.365  69.033  1.00 30.34           C  
ATOM    127  OE1 GLN A  94      37.589  22.914  69.386  1.00 31.20           O  
ATOM    128  NE2 GLN A  94      35.540  23.691  69.907  1.00 32.19           N  
ATOM    129  N   SER A  95      34.309  27.610  66.139  1.00 29.08           N  
ATOM    130  CA  SER A  95      33.295  28.660  66.259  1.00 28.73           C  
ATOM    131  C   SER A  95      32.714  29.081  64.920  1.00 28.63           C  
ATOM    132  O   SER A  95      31.659  29.719  64.861  1.00 28.34           O  
ATOM    133  CB  SER A  95      33.913  29.883  66.925  1.00 28.75           C  
ATOM    134  OG  SER A  95      34.869  30.490  66.068  1.00 29.17           O  
ATOM    135  N   LEU A  96      33.416  28.711  63.854  1.00 28.69           N  
ATOM    136  CA  LEU A  96      33.163  29.237  62.519  1.00 28.98           C  
ATOM    137  C   LEU A  96      31.902  28.715  61.871  1.00 28.83           C  
ATOM    138  O   LEU A  96      31.461  29.253  60.861  1.00 28.94           O  
ATOM    139  CB  LEU A  96      34.357  28.959  61.611  1.00 29.29           C  
ATOM    140  CG  LEU A  96      35.368  30.077  61.353  1.00 29.85           C  
ATOM    141  CD1 LEU A  96      35.520  31.041  62.522  1.00 29.91           C  
ATOM    142  CD2 LEU A  96      36.692  29.437  60.991  1.00 30.93           C  
ATOM    143  N   SER A  97      31.330  27.666  62.449  1.00 28.58           N  
ATOM    144  CA  SER A  97      30.036  27.163  62.008  1.00 28.55           C  
ATOM    145  C   SER A  97      28.895  27.984  62.623  1.00 28.14           C  
ATOM    146  O   SER A  97      27.717  27.719  62.362  1.00 27.91           O  
ATOM    147  CB  SER A  97      29.890  25.667  62.338  1.00 28.56           C  
ATOM    148  OG  SER A  97      30.042  25.423  63.729  1.00 29.32           O  
ATOM    149  N   GLY A  98      29.267  28.968  63.446  1.00 27.80           N  
ATOM    150  CA  GLY A  98      28.329  29.902  64.053  1.00 27.65           C  
ATOM    151  C   GLY A  98      28.202  31.157  63.212  1.00 27.93           C  
ATOM    152  O   GLY A  98      29.207  31.740  62.792  1.00 28.24           O  
ATOM    153  N   GLU A  99      26.960  31.577  62.976  1.00 27.88           N  
ATOM    154  CA  GLU A  99      26.649  32.721  62.118  1.00 27.92           C  
ATOM    155  C   GLU A  99      27.526  33.970  62.320  1.00 28.08           C  
ATOM    156  O   GLU A  99      28.258  34.375  61.410  1.00 28.33           O  
ATOM    157  CB  GLU A  99      25.167  33.077  62.248  1.00 27.85           C  
ATOM    158  CG  GLU A  99      24.717  34.189  61.327  1.00 27.43           C  
ATOM    159  CD  GLU A  99      23.218  34.229  61.148  1.00 28.66           C  
ATOM    160  OE1 GLU A  99      22.690  35.309  60.815  1.00 29.71           O  
ATOM    161  OE2 GLU A  99      22.557  33.184  61.328  1.00 29.52           O  
ATOM    162  N   SER A 100      27.438  34.575  63.502  1.00 28.00           N  
ATOM    163  CA  SER A 100      28.125  35.826  63.793  1.00 28.09           C  
ATOM    164  C   SER A 100      29.640  35.685  63.710  1.00 28.21           C  
ATOM    165  O   SER A 100      30.334  36.609  63.280  1.00 28.35           O  
ATOM    166  CB  SER A 100      27.731  36.337  65.177  1.00 28.05           C  
ATOM    167  OG  SER A 100      26.345  36.627  65.230  1.00 28.79           O  
ATOM    168  N   SER A 101      30.152  34.529  64.116  1.00 28.11           N  
ATOM    169  CA  SER A 101      31.592  34.312  64.110  1.00 28.12           C  
ATOM    170  C   SER A 101      32.088  34.180  62.681  1.00 28.07           C  
ATOM    171  O   SER A 101      33.150  34.710  62.340  1.00 28.22           O  
ATOM    172  CB  SER A 101      31.985  33.086  64.948  1.00 27.98           C  
ATOM    173  OG  SER A 101      31.818  33.338  66.336  1.00 27.56           O  
ATOM    174  N   ASN A 102      31.302  33.498  61.851  1.00 27.88           N  
ATOM    175  CA  ASN A 102      31.666  33.270  60.456  1.00 27.70           C  
ATOM    176  C   ASN A 102      31.709  34.583  59.683  1.00 27.89           C  
ATOM    177  O   ASN A 102      32.638  34.820  58.902  1.00 27.88           O  
ATOM    178  CB  ASN A 102      30.704  32.279  59.805  1.00 27.55           C  
ATOM    179  CG  ASN A 102      31.182  31.802  58.452  1.00 27.16           C  
ATOM    180  OD1 ASN A 102      31.003  32.479  57.446  1.00 26.64           O  
ATOM    181  ND2 ASN A 102      31.777  30.615  58.418  1.00 27.34           N  
ATOM    182  N   ARG A 103      30.720  35.442  59.928  1.00 27.77           N  
ATOM    183  CA  ARG A 103      30.630  36.733  59.242  1.00 27.79           C  
ATOM    184  C   ARG A 103      31.764  37.682  59.606  1.00 27.95           C  
ATOM    185  O   ARG A 103      32.302  38.363  58.742  1.00 28.07           O  
ATOM    186  CB  ARG A 103      29.281  37.405  59.509  1.00 27.46           C  
ATOM    187  CG  ARG A 103      28.122  36.605  58.996  1.00 27.21           C  
ATOM    188  CD  ARG A 103      27.080  37.489  58.385  1.00 25.25           C  
ATOM    189  NE  ARG A 103      26.019  36.685  57.800  1.00 23.85           N  
ATOM    190  CZ  ARG A 103      24.867  36.430  58.401  1.00 23.28           C  
ATOM    191  NH1 ARG A 103      24.635  36.931  59.605  1.00 24.56           N  
ATOM    192  NH2 ARG A 103      23.943  35.693  57.792  1.00 22.14           N  
ATOM    193  N   VAL A 104      32.110  37.733  60.887  1.00 28.17           N  
ATOM    194  CA  VAL A 104      33.182  38.596  61.359  1.00 28.75           C  
ATOM    195  C   VAL A 104      34.538  38.062  60.900  1.00 29.45           C  
ATOM    196  O   VAL A 104      35.439  38.832  60.585  1.00 29.36           O  
ATOM    197  CB  VAL A 104      33.116  38.794  62.911  1.00 28.63           C  
ATOM    198  CG1 VAL A 104      34.477  39.132  63.508  1.00 27.61           C  
ATOM    199  CG2 VAL A 104      32.085  39.869  63.264  1.00 28.43           C  
HETATM  200  N   MSE A 105      34.669  36.740  60.846  1.00 30.64           N  
HETATM  201  CA  MSE A 105      35.918  36.123  60.417  1.00 31.50           C  
HETATM  202  C   MSE A 105      36.241  36.413  58.950  1.00 30.93           C  
HETATM  203  O   MSE A 105      37.374  36.770  58.637  1.00 30.87           O  
HETATM  204  CB  MSE A 105      35.912  34.610  60.667  1.00 32.75           C  
HETATM  205  CG  MSE A 105      37.295  33.974  60.494  1.00 35.49           C  
HETATM  206 SE   MSE A 105      38.579  34.624  61.820  1.00 45.42          SE  
HETATM  207  CE  MSE A 105      40.000  35.200  60.616  1.00 38.71           C  
ATOM    208  N   TRP A 106      35.240  36.278  58.075  1.00 30.06           N  
ATOM    209  CA  TRP A 106      35.467  36.221  56.632  1.00 29.43           C  
ATOM    210  C   TRP A 106      34.893  37.404  55.837  1.00 29.67           C  
ATOM    211  O   TRP A 106      34.721  37.320  54.611  1.00 29.78           O  
ATOM    212  CB  TRP A 106      34.963  34.887  56.044  1.00 28.62           C  
ATOM    213  CG  TRP A 106      35.456  33.624  56.739  1.00 27.93           C  
ATOM    214  CD1 TRP A 106      34.684  32.587  57.187  1.00 26.75           C  
ATOM    215  CD2 TRP A 106      36.818  33.269  57.048  1.00 27.17           C  
ATOM    216  NE1 TRP A 106      35.471  31.618  57.754  1.00 26.78           N  
ATOM    217  CE2 TRP A 106      36.784  32.006  57.682  1.00 27.10           C  
ATOM    218  CE3 TRP A 106      38.061  33.897  56.853  1.00 26.49           C  
ATOM    219  CZ2 TRP A 106      37.948  31.356  58.128  1.00 26.97           C  
ATOM    220  CZ3 TRP A 106      39.220  33.250  57.300  1.00 27.13           C  
ATOM    221  CH2 TRP A 106      39.150  31.993  57.930  1.00 27.20           C  
ATOM    222  N   ASN A 107      34.599  38.502  56.521  1.00 29.84           N  
ATOM    223  CA  ASN A 107      34.313  39.762  55.836  1.00 29.92           C  
ATOM    224  C   ASN A 107      35.556  40.169  55.034  1.00 30.44           C  
ATOM    225  O   ASN A 107      36.684  39.863  55.443  1.00 30.55           O  
ATOM    226  CB  ASN A 107      33.954  40.838  56.857  1.00 29.41           C  
ATOM    227  CG  ASN A 107      33.306  42.053  56.229  1.00 29.18           C  
ATOM    228  OD1 ASN A 107      33.972  43.059  55.947  1.00 28.47           O  
ATOM    229  ND2 ASN A 107      31.999  41.978  56.018  1.00 27.76           N  
ATOM    230  N   ASP A 108      35.368  40.833  53.895  1.00 30.88           N  
ATOM    231  CA  ASP A 108      36.518  41.240  53.070  1.00 31.48           C  
ATOM    232  C   ASP A 108      36.856  42.748  53.102  1.00 31.63           C  
ATOM    233  O   ASP A 108      37.597  43.247  52.244  1.00 31.49           O  
ATOM    234  CB  ASP A 108      36.356  40.745  51.630  1.00 31.75           C  
ATOM    235  CG  ASP A 108      35.130  41.316  50.947  1.00 32.90           C  
ATOM    236  OD1 ASP A 108      34.272  41.914  51.631  1.00 34.90           O  
ATOM    237  OD2 ASP A 108      35.020  41.159  49.715  1.00 34.75           O  
ATOM    238  N   ARG A 109      36.313  43.466  54.084  1.00 31.71           N  
ATOM    239  CA  ARG A 109      36.756  44.830  54.344  1.00 31.85           C  
ATOM    240  C   ARG A 109      37.994  44.768  55.217  1.00 31.38           C  
ATOM    241  O   ARG A 109      37.955  44.247  56.334  1.00 31.47           O  
ATOM    242  CB  ARG A 109      35.671  45.666  55.023  1.00 32.30           C  
ATOM    243  CG  ARG A 109      36.138  47.077  55.353  1.00 33.49           C  
ATOM    244  CD  ARG A 109      35.230  47.764  56.356  1.00 35.64           C  
ATOM    245  NE  ARG A 109      34.214  48.597  55.709  1.00 36.26           N  
ATOM    246  CZ  ARG A 109      34.420  49.834  55.260  1.00 36.05           C  
ATOM    247  NH1 ARG A 109      35.613  50.414  55.363  1.00 35.49           N  
ATOM    248  NH2 ARG A 109      33.420  50.492  54.701  1.00 36.99           N  
ATOM    249  N   TYR A 110      39.086  45.315  54.704  1.00 30.93           N  
ATOM    250  CA  TYR A 110      40.390  45.146  55.323  1.00 30.70           C  
ATOM    251  C   TYR A 110      40.940  46.416  55.939  1.00 30.61           C  
ATOM    252  O   TYR A 110      42.070  46.430  56.423  1.00 31.03           O  
ATOM    253  CB  TYR A 110      41.394  44.626  54.293  1.00 30.38           C  
ATOM    254  CG  TYR A 110      40.998  43.320  53.660  1.00 30.24           C  
ATOM    255  CD1 TYR A 110      40.649  42.223  54.447  1.00 29.49           C  
ATOM    256  CD2 TYR A 110      40.988  43.171  52.274  1.00 29.78           C  
ATOM    257  CE1 TYR A 110      40.288  41.017  53.875  1.00 29.60           C  
ATOM    258  CE2 TYR A 110      40.626  41.964  51.690  1.00 29.86           C  
ATOM    259  CZ  TYR A 110      40.279  40.895  52.504  1.00 29.92           C  
ATOM    260  OH  TYR A 110      39.926  39.696  51.952  1.00 30.80           O  
ATOM    261  N   ASP A 111      40.165  47.487  55.913  1.00 30.31           N  
ATOM    262  CA  ASP A 111      40.730  48.779  56.288  1.00 30.08           C  
ATOM    263  C   ASP A 111      40.358  49.215  57.714  1.00 29.81           C  
ATOM    264  O   ASP A 111      40.858  50.228  58.210  1.00 29.97           O  
ATOM    265  CB  ASP A 111      40.433  49.855  55.219  1.00 30.12           C  
ATOM    266  CG  ASP A 111      38.951  50.177  55.084  1.00 29.30           C  
ATOM    267  OD1 ASP A 111      38.160  49.299  54.700  1.00 28.08           O  
ATOM    268  OD2 ASP A 111      38.583  51.340  55.341  1.00 29.97           O  
ATOM    269  N   THR A 112      39.497  48.430  58.363  1.00 29.24           N  
ATOM    270  CA  THR A 112      39.196  48.592  59.786  1.00 28.46           C  
ATOM    271  C   THR A 112      39.624  47.362  60.603  1.00 27.82           C  
ATOM    272  O   THR A 112      39.813  46.276  60.062  1.00 27.21           O  
ATOM    273  CB  THR A 112      37.694  48.913  60.046  1.00 28.64           C  
ATOM    274  OG1 THR A 112      36.946  47.702  60.281  1.00 27.97           O  
ATOM    275  CG2 THR A 112      37.106  49.708  58.886  1.00 28.60           C  
ATOM    276  N   LEU A 113      39.775  47.567  61.908  1.00 27.37           N  
ATOM    277  CA  LEU A 113      40.154  46.529  62.853  1.00 27.08           C  
ATOM    278  C   LEU A 113      39.178  45.373  62.711  1.00 27.18           C  
ATOM    279  O   LEU A 113      37.963  45.588  62.806  1.00 27.78           O  
ATOM    280  CB  LEU A 113      40.082  47.104  64.269  1.00 26.74           C  
ATOM    281  CG  LEU A 113      41.191  46.951  65.314  1.00 26.74           C  
ATOM    282  CD1 LEU A 113      42.611  46.825  64.745  1.00 26.52           C  
ATOM    283  CD2 LEU A 113      41.105  48.125  66.259  1.00 25.90           C  
ATOM    284  N   LEU A 114      39.702  44.168  62.455  1.00 26.55           N  
ATOM    285  CA  LEU A 114      38.879  42.967  62.310  1.00 26.11           C  
ATOM    286  C   LEU A 114      37.905  42.817  63.481  1.00 25.96           C  
ATOM    287  O   LEU A 114      36.710  42.620  63.269  1.00 25.78           O  
ATOM    288  CB  LEU A 114      39.745  41.700  62.151  1.00 26.07           C  
ATOM    289  CG  LEU A 114      39.027  40.352  61.880  1.00 26.54           C  
ATOM    290  CD1 LEU A 114      39.908  39.367  61.104  1.00 25.68           C  
ATOM    291  CD2 LEU A 114      38.483  39.687  63.151  1.00 25.72           C  
ATOM    292  N   ILE A 115      38.424  42.916  64.709  1.00 25.71           N  
ATOM    293  CA  ILE A 115      37.609  42.796  65.929  1.00 25.22           C  
ATOM    294  C   ILE A 115      36.576  43.924  66.073  1.00 24.95           C  
ATOM    295  O   ILE A 115      35.705  43.872  66.951  1.00 24.85           O  
ATOM    296  CB  ILE A 115      38.483  42.693  67.220  1.00 25.16           C  
ATOM    297  CG1 ILE A 115      39.265  43.992  67.480  1.00 25.76           C  
ATOM    298  CG2 ILE A 115      39.424  41.479  67.158  1.00 25.06           C  
ATOM    299  CD1 ILE A 115      39.788  44.138  68.938  1.00 25.61           C  
ATOM    300  N   ALA A 116      36.681  44.943  65.219  1.00 24.54           N  
ATOM    301  CA  ALA A 116      35.711  46.039  65.212  1.00 24.35           C  
ATOM    302  C   ALA A 116      34.504  45.745  64.303  1.00 24.12           C  
ATOM    303  O   ALA A 116      33.514  46.484  64.323  1.00 24.25           O  
ATOM    304  CB  ALA A 116      36.379  47.359  64.831  1.00 23.99           C  
ATOM    305  N   ARG A 117      34.598  44.673  63.517  1.00 23.57           N  
ATOM    306  CA  ARG A 117      33.512  44.234  62.657  1.00 23.31           C  
ATOM    307  C   ARG A 117      32.309  43.847  63.508  1.00 23.56           C  
ATOM    308  O   ARG A 117      32.425  43.049  64.439  1.00 23.46           O  
ATOM    309  CB  ARG A 117      33.941  43.042  61.800  1.00 23.21           C  
ATOM    310  CG  ARG A 117      34.941  43.355  60.696  1.00 22.37           C  
ATOM    311  CD  ARG A 117      35.543  42.069  60.146  1.00 21.28           C  
ATOM    312  NE  ARG A 117      36.586  42.321  59.154  1.00 21.17           N  
ATOM    313  CZ  ARG A 117      37.237  41.379  58.468  1.00 20.82           C  
ATOM    314  NH1 ARG A 117      36.967  40.088  58.648  1.00 19.83           N  
ATOM    315  NH2 ARG A 117      38.167  41.734  57.587  1.00 19.97           N  
ATOM    316  N   ASP A 118      31.160  44.430  63.184  1.00 23.74           N  
ATOM    317  CA  ASP A 118      29.932  44.201  63.924  1.00 23.87           C  
ATOM    318  C   ASP A 118      29.037  43.229  63.159  1.00 24.17           C  
ATOM    319  O   ASP A 118      28.662  43.509  62.021  1.00 24.12           O  
ATOM    320  CB  ASP A 118      29.210  45.535  64.147  1.00 23.84           C  
ATOM    321  CG  ASP A 118      27.910  45.380  64.915  1.00 23.54           C  
ATOM    322  OD1 ASP A 118      27.832  44.552  65.835  1.00 22.31           O  
ATOM    323  OD2 ASP A 118      26.956  46.101  64.602  1.00 25.28           O  
ATOM    324  N   PRO A 119      28.673  42.095  63.792  1.00 24.60           N  
ATOM    325  CA  PRO A 119      27.903  41.019  63.150  1.00 25.01           C  
ATOM    326  C   PRO A 119      26.496  41.439  62.707  1.00 25.46           C  
ATOM    327  O   PRO A 119      25.959  40.864  61.762  1.00 25.50           O  
ATOM    328  CB  PRO A 119      27.802  39.951  64.242  1.00 25.25           C  
ATOM    329  CG  PRO A 119      28.832  40.315  65.257  1.00 25.12           C  
ATOM    330  CD  PRO A 119      28.971  41.785  65.201  1.00 24.54           C  
ATOM    331  N   ARG A 120      25.903  42.416  63.385  1.00 25.84           N  
ATOM    332  CA  ARG A 120      24.624  42.960  62.937  1.00 26.53           C  
ATOM    333  C   ARG A 120      24.787  44.027  61.852  1.00 26.57           C  
ATOM    334  O   ARG A 120      23.924  44.137  60.981  1.00 26.52           O  
ATOM    335  CB  ARG A 120      23.761  43.466  64.103  1.00 26.45           C  
ATOM    336  CG  ARG A 120      24.535  43.889  65.333  1.00 28.54           C  
ATOM    337  CD  ARG A 120      23.987  45.182  65.936  1.00 31.82           C  
ATOM    338  NE  ARG A 120      23.867  46.236  64.920  1.00 33.39           N  
ATOM    339  CZ  ARG A 120      24.065  47.537  65.128  1.00 33.91           C  
ATOM    340  NH1 ARG A 120      23.923  48.382  64.118  1.00 34.52           N  
ATOM    341  NH2 ARG A 120      24.414  48.000  66.327  1.00 33.76           N  
ATOM    342  N   GLU A 121      25.878  44.801  61.898  1.00 26.72           N  
ATOM    343  CA  GLU A 121      26.178  45.776  60.834  1.00 27.37           C  
ATOM    344  C   GLU A 121      26.365  45.071  59.499  1.00 27.28           C  
ATOM    345  O   GLU A 121      25.807  45.495  58.485  1.00 27.28           O  
ATOM    346  CB  GLU A 121      27.426  46.614  61.144  1.00 27.72           C  
ATOM    347  CG  GLU A 121      27.721  47.729  60.100  1.00 30.02           C  
ATOM    348  CD  GLU A 121      28.940  47.466  59.178  1.00 32.64           C  
ATOM    349  OE1 GLU A 121      29.886  48.292  59.218  1.00 33.31           O  
ATOM    350  OE2 GLU A 121      28.957  46.470  58.407  1.00 32.54           O  
ATOM    351  N   ILE A 122      27.156  43.997  59.515  1.00 27.19           N  
ATOM    352  CA  ILE A 122      27.368  43.163  58.342  1.00 27.06           C  
ATOM    353  C   ILE A 122      26.058  42.526  57.871  1.00 27.39           C  
ATOM    354  O   ILE A 122      25.789  42.493  56.681  1.00 27.33           O  
ATOM    355  CB  ILE A 122      28.474  42.109  58.579  1.00 26.77           C  
ATOM    356  CG1 ILE A 122      29.792  42.808  58.925  1.00 26.46           C  
ATOM    357  CG2 ILE A 122      28.667  41.253  57.338  1.00 26.40           C  
ATOM    358  CD1 ILE A 122      30.787  41.963  59.727  1.00 24.92           C  
ATOM    359  N   LYS A 123      25.235  42.045  58.797  1.00 27.98           N  
ATOM    360  CA  LYS A 123      23.956  41.435  58.426  1.00 28.55           C  
ATOM    361  C   LYS A 123      23.006  42.448  57.795  1.00 28.76           C  
ATOM    362  O   LYS A 123      22.312  42.133  56.819  1.00 28.98           O  
ATOM    363  CB  LYS A 123      23.301  40.723  59.620  1.00 29.04           C  
ATOM    364  CG  LYS A 123      21.895  40.182  59.347  1.00 30.13           C  
ATOM    365  CD  LYS A 123      21.721  38.763  59.836  1.00 31.98           C  
ATOM    366  CE  LYS A 123      20.251  38.341  59.800  1.00 33.31           C  
ATOM    367  NZ  LYS A 123      20.056  37.017  60.475  1.00 33.26           N  
ATOM    368  N   ASN A 124      22.984  43.664  58.339  1.00 28.92           N  
ATOM    369  CA  ASN A 124      22.187  44.734  57.748  1.00 29.12           C  
ATOM    370  C   ASN A 124      22.662  45.055  56.348  1.00 29.00           C  
ATOM    371  O   ASN A 124      21.848  45.218  55.442  1.00 28.90           O  
ATOM    372  CB  ASN A 124      22.176  45.983  58.626  1.00 29.26           C  
ATOM    373  CG  ASN A 124      21.166  45.884  59.769  1.00 30.45           C  
ATOM    374  OD1 ASN A 124      20.704  44.791  60.126  1.00 31.06           O  
ATOM    375  ND2 ASN A 124      20.820  47.031  60.348  1.00 31.53           N  
ATOM    376  N   ALA A 125      23.982  45.112  56.177  1.00 29.08           N  
ATOM    377  CA  ALA A 125      24.597  45.345  54.872  1.00 29.23           C  
ATOM    378  C   ALA A 125      24.147  44.312  53.831  1.00 29.49           C  
ATOM    379  O   ALA A 125      23.803  44.679  52.703  1.00 29.94           O  
ATOM    380  CB  ALA A 125      26.118  45.399  54.985  1.00 28.81           C  
ATOM    381  N   ILE A 126      24.117  43.035  54.210  1.00 29.54           N  
ATOM    382  CA  ILE A 126      23.615  41.990  53.315  1.00 29.83           C  
ATOM    383  C   ILE A 126      22.122  42.202  52.989  1.00 30.19           C  
ATOM    384  O   ILE A 126      21.733  42.168  51.812  1.00 30.18           O  
ATOM    385  CB  ILE A 126      23.855  40.551  53.858  1.00 29.71           C  
ATOM    386  CG1 ILE A 126      25.342  40.317  54.161  1.00 28.82           C  
ATOM    387  CG2 ILE A 126      23.347  39.513  52.852  1.00 29.42           C  
ATOM    388  CD1 ILE A 126      25.620  39.234  55.187  1.00 26.42           C  
ATOM    389  N   GLU A 127      21.304  42.429  54.018  1.00 30.29           N  
ATOM    390  CA  GLU A 127      19.874  42.699  53.814  1.00 31.14           C  
ATOM    391  C   GLU A 127      19.669  43.837  52.816  1.00 30.91           C  
ATOM    392  O   GLU A 127      18.807  43.754  51.935  1.00 30.86           O  
ATOM    393  CB  GLU A 127      19.179  43.045  55.131  1.00 30.87           C  
ATOM    394  CG  GLU A 127      18.974  41.860  56.058  1.00 32.06           C  
ATOM    395  CD  GLU A 127      18.642  42.274  57.494  1.00 32.74           C  
ATOM    396  OE1 GLU A 127      18.416  43.483  57.755  1.00 34.51           O  
ATOM    397  OE2 GLU A 127      18.604  41.377  58.365  1.00 34.93           O  
ATOM    398  N   LYS A 128      20.481  44.885  52.962  1.00 30.77           N  
ATOM    399  CA  LYS A 128      20.406  46.058  52.108  1.00 30.76           C  
ATOM    400  C   LYS A 128      20.754  45.654  50.684  1.00 30.66           C  
ATOM    401  O   LYS A 128      19.952  45.837  49.759  1.00 30.91           O  
ATOM    402  CB  LYS A 128      21.345  47.156  52.619  1.00 30.62           C  
ATOM    403  CG  LYS A 128      21.368  48.445  51.798  1.00 30.91           C  
ATOM    404  CD  LYS A 128      19.996  49.142  51.730  1.00 31.59           C  
ATOM    405  CE  LYS A 128      20.125  50.631  51.394  1.00 29.99           C  
ATOM    406  NZ  LYS A 128      20.832  50.864  50.109  1.00 29.77           N  
ATOM    407  N   SER A 129      21.940  45.076  50.530  1.00 30.08           N  
ATOM    408  CA  SER A 129      22.414  44.602  49.249  1.00 29.89           C  
ATOM    409  C   SER A 129      21.348  43.760  48.526  1.00 29.52           C  
ATOM    410  O   SER A 129      21.063  43.988  47.353  1.00 29.28           O  
ATOM    411  CB  SER A 129      23.712  43.817  49.444  1.00 29.68           C  
ATOM    412  OG  SER A 129      24.046  43.119  48.257  1.00 31.02           O  
ATOM    413  N   VAL A 130      20.755  42.809  49.242  1.00 29.34           N  
ATOM    414  CA  VAL A 130      19.672  41.978  48.705  1.00 29.10           C  
ATOM    415  C   VAL A 130      18.448  42.821  48.314  1.00 29.39           C  
ATOM    416  O   VAL A 130      17.926  42.668  47.198  1.00 29.77           O  
ATOM    417  CB  VAL A 130      19.288  40.822  49.682  1.00 29.04           C  
ATOM    418  CG1 VAL A 130      17.971  40.198  49.307  1.00 28.21           C  
ATOM    419  CG2 VAL A 130      20.374  39.759  49.703  1.00 28.49           C  
ATOM    420  N   THR A 131      18.011  43.722  49.200  1.00 28.89           N  
ATOM    421  CA  THR A 131      16.904  44.630  48.875  1.00 28.66           C  
ATOM    422  C   THR A 131      17.182  45.379  47.559  1.00 28.56           C  
ATOM    423  O   THR A 131      16.315  45.440  46.681  1.00 28.23           O  
ATOM    424  CB  THR A 131      16.584  45.636  50.035  1.00 28.75           C  
ATOM    425  OG1 THR A 131      16.291  44.916  51.237  1.00 28.64           O  
ATOM    426  CG2 THR A 131      15.386  46.524  49.692  1.00 28.08           C  
ATOM    427  N   ASP A 132      18.395  45.921  47.429  1.00 28.56           N  
ATOM    428  CA  ASP A 132      18.816  46.645  46.218  1.00 28.76           C  
ATOM    429  C   ASP A 132      18.645  45.809  44.949  1.00 28.81           C  
ATOM    430  O   ASP A 132      18.329  46.349  43.886  1.00 28.78           O  
ATOM    431  CB  ASP A 132      20.275  47.107  46.330  1.00 28.74           C  
ATOM    432  CG  ASP A 132      20.459  48.288  47.276  1.00 29.51           C  
ATOM    433  OD1 ASP A 132      19.453  48.894  47.712  1.00 29.94           O  
ATOM    434  OD2 ASP A 132      21.630  48.616  47.585  1.00 31.01           O  
ATOM    435  N   PHE A 133      18.854  44.494  45.068  1.00 28.69           N  
ATOM    436  CA  PHE A 133      18.703  43.570  43.934  1.00 28.41           C  
ATOM    437  C   PHE A 133      17.280  43.006  43.791  1.00 28.11           C  
ATOM    438  O   PHE A 133      17.063  42.056  43.039  1.00 27.86           O  
ATOM    439  CB  PHE A 133      19.702  42.418  44.043  1.00 28.20           C  
ATOM    440  CG  PHE A 133      21.072  42.740  43.531  1.00 28.45           C  
ATOM    441  CD1 PHE A 133      21.394  42.539  42.192  1.00 28.59           C  
ATOM    442  CD2 PHE A 133      22.060  43.209  44.387  1.00 28.70           C  
ATOM    443  CE1 PHE A 133      22.678  42.813  41.711  1.00 27.92           C  
ATOM    444  CE2 PHE A 133      23.353  43.479  43.906  1.00 28.41           C  
ATOM    445  CZ  PHE A 133      23.655  43.286  42.569  1.00 27.01           C  
ATOM    446  N   GLY A 134      16.326  43.592  44.512  1.00 27.84           N  
ATOM    447  CA  GLY A 134      14.916  43.215  44.401  1.00 27.57           C  
ATOM    448  C   GLY A 134      14.557  41.883  45.033  1.00 27.47           C  
ATOM    449  O   GLY A 134      13.625  41.220  44.591  1.00 27.36           O  
ATOM    450  N   GLY A 135      15.296  41.491  46.070  1.00 27.62           N  
ATOM    451  CA  GLY A 135      14.999  40.277  46.826  1.00 27.68           C  
ATOM    452  C   GLY A 135      16.023  39.179  46.644  1.00 28.12           C  
ATOM    453  O   GLY A 135      16.789  39.181  45.674  1.00 28.06           O  
ATOM    454  N   LEU A 136      16.018  38.231  47.580  1.00 28.71           N  
ATOM    455  CA  LEU A 136      16.946  37.091  47.586  1.00 29.42           C  
ATOM    456  C   LEU A 136      16.899  36.219  46.325  1.00 29.81           C  
ATOM    457  O   LEU A 136      17.946  35.880  45.762  1.00 29.99           O  
ATOM    458  CB  LEU A 136      16.721  36.229  48.834  1.00 29.15           C  
ATOM    459  CG  LEU A 136      17.740  35.134  49.171  1.00 29.23           C  
ATOM    460  CD1 LEU A 136      19.168  35.660  49.206  1.00 29.50           C  
ATOM    461  CD2 LEU A 136      17.393  34.502  50.508  1.00 29.62           C  
ATOM    462  N   GLU A 137      15.698  35.849  45.892  1.00 30.27           N  
ATOM    463  CA  GLU A 137      15.543  35.077  44.665  1.00 31.12           C  
ATOM    464  C   GLU A 137      16.313  35.758  43.529  1.00 31.10           C  
ATOM    465  O   GLU A 137      17.158  35.137  42.878  1.00 31.03           O  
ATOM    466  CB  GLU A 137      14.066  34.959  44.284  1.00 31.52           C  
ATOM    467  CG  GLU A 137      13.166  34.317  45.330  1.00 33.24           C  
ATOM    468  CD  GLU A 137      13.024  32.822  45.127  1.00 35.44           C  
ATOM    469  OE1 GLU A 137      11.892  32.358  44.803  1.00 35.50           O  
ATOM    470  OE2 GLU A 137      14.055  32.123  45.282  1.00 35.79           O  
ATOM    471  N   ASN A 138      16.019  37.041  43.315  1.00 31.15           N  
ATOM    472  CA  ASN A 138      16.632  37.817  42.244  1.00 31.18           C  
ATOM    473  C   ASN A 138      18.120  38.054  42.442  1.00 31.07           C  
ATOM    474  O   ASN A 138      18.882  38.102  41.470  1.00 30.75           O  
ATOM    475  CB  ASN A 138      15.904  39.145  42.061  1.00 31.20           C  
ATOM    476  CG  ASN A 138      15.268  39.265  40.706  1.00 31.93           C  
ATOM    477  OD1 ASN A 138      14.087  38.959  40.532  1.00 32.84           O  
ATOM    478  ND2 ASN A 138      16.055  39.691  39.719  1.00 32.56           N  
ATOM    479  N   TYR A 139      18.525  38.212  43.698  1.00 30.97           N  
ATOM    480  CA  TYR A 139      19.928  38.364  44.005  1.00 31.38           C  
ATOM    481  C   TYR A 139      20.652  37.102  43.551  1.00 31.79           C  
ATOM    482  O   TYR A 139      21.627  37.181  42.817  1.00 31.70           O  
ATOM    483  CB  TYR A 139      20.161  38.620  45.500  1.00 31.30           C  
ATOM    484  CG  TYR A 139      21.625  38.797  45.854  1.00 31.14           C  
ATOM    485  CD1 TYR A 139      22.179  40.067  45.974  1.00 31.36           C  
ATOM    486  CD2 TYR A 139      22.461  37.693  46.052  1.00 30.75           C  
ATOM    487  CE1 TYR A 139      23.527  40.242  46.294  1.00 31.72           C  
ATOM    488  CE2 TYR A 139      23.808  37.854  46.367  1.00 31.07           C  
ATOM    489  CZ  TYR A 139      24.334  39.133  46.486  1.00 31.72           C  
ATOM    490  OH  TYR A 139      25.663  39.309  46.798  1.00 32.04           O  
ATOM    491  N   LYS A 140      20.151  35.945  43.972  1.00 32.25           N  
ATOM    492  CA  LYS A 140      20.757  34.673  43.622  1.00 32.92           C  
ATOM    493  C   LYS A 140      20.857  34.499  42.106  1.00 33.36           C  
ATOM    494  O   LYS A 140      21.888  34.063  41.586  1.00 33.32           O  
ATOM    495  CB  LYS A 140      19.977  33.521  44.263  1.00 33.00           C  
ATOM    496  CG  LYS A 140      20.197  33.407  45.768  1.00 33.27           C  
ATOM    497  CD  LYS A 140      19.085  32.636  46.443  1.00 33.48           C  
ATOM    498  CE  LYS A 140      19.499  31.222  46.735  1.00 33.94           C  
ATOM    499  NZ  LYS A 140      18.351  30.426  47.236  1.00 35.20           N  
ATOM    500  N   GLU A 141      19.793  34.872  41.401  1.00 33.93           N  
ATOM    501  CA  GLU A 141      19.734  34.717  39.953  1.00 34.44           C  
ATOM    502  C   GLU A 141      20.794  35.555  39.234  1.00 34.64           C  
ATOM    503  O   GLU A 141      21.417  35.090  38.278  1.00 34.52           O  
ATOM    504  CB  GLU A 141      18.338  35.063  39.453  1.00 34.52           C  
ATOM    505  CG  GLU A 141      17.959  34.379  38.154  1.00 35.31           C  
ATOM    506  CD  GLU A 141      16.468  34.136  38.058  1.00 35.70           C  
ATOM    507  OE1 GLU A 141      15.972  33.225  38.759  1.00 35.87           O  
ATOM    508  OE2 GLU A 141      15.797  34.855  37.287  1.00 35.22           O  
ATOM    509  N   LEU A 142      21.003  36.778  39.717  1.00 34.94           N  
ATOM    510  CA  LEU A 142      21.965  37.710  39.120  1.00 35.36           C  
ATOM    511  C   LEU A 142      23.421  37.500  39.579  1.00 35.55           C  
ATOM    512  O   LEU A 142      24.352  37.898  38.876  1.00 35.79           O  
ATOM    513  CB  LEU A 142      21.528  39.158  39.377  1.00 35.14           C  
ATOM    514  CG  LEU A 142      20.576  39.908  38.421  1.00 35.62           C  
ATOM    515  CD1 LEU A 142      19.291  39.159  38.058  1.00 36.01           C  
ATOM    516  CD2 LEU A 142      20.213  41.260  39.016  1.00 35.58           C  
ATOM    517  N   THR A 143      23.610  36.858  40.734  1.00 35.69           N  
ATOM    518  CA  THR A 143      24.923  36.785  41.390  1.00 35.80           C  
ATOM    519  C   THR A 143      25.608  35.414  41.323  1.00 35.53           C  
ATOM    520  O   THR A 143      26.698  35.233  41.884  1.00 35.62           O  
ATOM    521  CB  THR A 143      24.839  37.210  42.875  1.00 36.00           C  
ATOM    522  OG1 THR A 143      23.686  38.032  43.087  1.00 36.82           O  
ATOM    523  CG2 THR A 143      26.070  38.001  43.276  1.00 36.98           C  
ATOM    524  N   GLY A 144      24.984  34.457  40.639  1.00 35.26           N  
ATOM    525  CA  GLY A 144      25.561  33.119  40.488  1.00 34.68           C  
ATOM    526  C   GLY A 144      24.819  32.010  41.215  1.00 34.49           C  
ATOM    527  O   GLY A 144      25.007  30.827  40.911  1.00 34.52           O  
ATOM    528  N   GLY A 145      23.983  32.393  42.180  1.00 34.24           N  
ATOM    529  CA  GLY A 145      23.145  31.457  42.926  1.00 33.64           C  
ATOM    530  C   GLY A 145      23.421  31.486  44.417  1.00 33.47           C  
ATOM    531  O   GLY A 145      22.730  30.824  45.196  1.00 33.07           O  
ATOM    532  N   ALA A 146      24.430  32.260  44.813  1.00 33.46           N  
ATOM    533  CA  ALA A 146      24.881  32.289  46.207  1.00 33.45           C  
ATOM    534  C   ALA A 146      23.883  32.983  47.115  1.00 33.12           C  
ATOM    535  O   ALA A 146      23.421  34.076  46.811  1.00 33.08           O  
ATOM    536  CB  ALA A 146      26.256  32.956  46.322  1.00 33.66           C  
ATOM    537  N   ASP A 147      23.537  32.323  48.217  1.00 32.80           N  
ATOM    538  CA  ASP A 147      22.823  32.965  49.310  1.00 32.39           C  
ATOM    539  C   ASP A 147      23.870  33.686  50.152  1.00 32.07           C  
ATOM    540  O   ASP A 147      24.696  33.037  50.798  1.00 31.94           O  
ATOM    541  CB  ASP A 147      22.084  31.926  50.153  1.00 32.42           C  
ATOM    542  CG  ASP A 147      21.133  32.552  51.169  1.00 33.12           C  
ATOM    543  OD1 ASP A 147      21.498  33.565  51.808  1.00 32.82           O  
ATOM    544  OD2 ASP A 147      20.010  32.017  51.331  1.00 34.47           O  
ATOM    545  N   PRO A 148      23.838  35.030  50.156  1.00 31.97           N  
ATOM    546  CA  PRO A 148      24.883  35.798  50.838  1.00 31.89           C  
ATOM    547  C   PRO A 148      24.725  35.798  52.366  1.00 32.02           C  
ATOM    548  O   PRO A 148      25.550  36.390  53.070  1.00 32.32           O  
ATOM    549  CB  PRO A 148      24.711  37.208  50.257  1.00 31.80           C  
ATOM    550  CG  PRO A 148      23.236  37.307  49.957  1.00 32.19           C  
ATOM    551  CD  PRO A 148      22.809  35.904  49.551  1.00 32.23           C  
ATOM    552  N   PHE A 149      23.679  35.135  52.861  1.00 31.91           N  
ATOM    553  CA  PHE A 149      23.466  34.957  54.296  1.00 31.86           C  
ATOM    554  C   PHE A 149      24.151  33.707  54.835  1.00 31.55           C  
ATOM    555  O   PHE A 149      24.309  33.557  56.043  1.00 31.49           O  
ATOM    556  CB  PHE A 149      21.970  34.874  54.618  1.00 31.99           C  
ATOM    557  CG  PHE A 149      21.222  36.167  54.416  1.00 32.37           C  
ATOM    558  CD1 PHE A 149      21.381  37.229  55.303  1.00 32.45           C  
ATOM    559  CD2 PHE A 149      20.337  36.313  53.352  1.00 32.62           C  
ATOM    560  CE1 PHE A 149      20.682  38.418  55.121  1.00 32.20           C  
ATOM    561  CE2 PHE A 149      19.631  37.502  53.165  1.00 32.41           C  
ATOM    562  CZ  PHE A 149      19.808  38.554  54.050  1.00 31.88           C  
ATOM    563  N   ALA A 150      24.545  32.811  53.938  1.00 31.44           N  
ATOM    564  CA  ALA A 150      25.086  31.509  54.319  1.00 31.37           C  
ATOM    565  C   ALA A 150      26.496  31.596  54.902  1.00 31.51           C  
ATOM    566  O   ALA A 150      27.240  32.530  54.610  1.00 31.55           O  
ATOM    567  CB  ALA A 150      25.066  30.571  53.129  1.00 31.18           C  
ATOM    568  N   LEU A 151      26.841  30.619  55.736  1.00 31.65           N  
ATOM    569  CA  LEU A 151      28.199  30.463  56.248  1.00 32.10           C  
ATOM    570  C   LEU A 151      29.180  30.295  55.096  1.00 32.11           C  
ATOM    571  O   LEU A 151      28.922  29.542  54.160  1.00 31.94           O  
ATOM    572  CB  LEU A 151      28.286  29.241  57.175  1.00 32.07           C  
ATOM    573  CG  LEU A 151      27.315  29.124  58.361  1.00 32.42           C  
ATOM    574  CD1 LEU A 151      27.334  27.697  58.941  1.00 31.91           C  
ATOM    575  CD2 LEU A 151      27.596  30.171  59.443  1.00 31.74           C  
HETATM  576  N   MSE A 152      30.306  30.990  55.172  1.00 32.61           N  
HETATM  577  CA  MSE A 152      31.283  30.975  54.080  1.00 33.27           C  
HETATM  578  C   MSE A 152      32.404  29.942  54.246  1.00 32.14           C  
HETATM  579  O   MSE A 152      32.920  29.743  55.341  1.00 31.69           O  
HETATM  580  CB  MSE A 152      31.872  32.378  53.845  1.00 34.81           C  
HETATM  581  CG  MSE A 152      30.847  33.428  53.372  1.00 38.12           C  
HETATM  582 SE   MSE A 152      29.635  32.715  51.982  1.00 50.22          SE  
HETATM  583  CE  MSE A 152      28.188  34.056  52.014  1.00 42.03           C  
ATOM    584  N   THR A 153      32.733  29.275  53.141  1.00 31.06           N  
ATOM    585  CA  THR A 153      33.973  28.511  52.989  1.00 29.78           C  
ATOM    586  C   THR A 153      34.947  29.423  52.231  1.00 28.65           C  
ATOM    587  O   THR A 153      34.843  29.574  51.010  1.00 28.40           O  
ATOM    588  CB  THR A 153      33.748  27.183  52.198  1.00 29.86           C  
ATOM    589  OG1 THR A 153      32.749  26.386  52.849  1.00 30.16           O  
ATOM    590  CG2 THR A 153      35.041  26.364  52.084  1.00 29.45           C  
ATOM    591  N   PRO A 154      35.887  30.046  52.960  1.00 27.78           N  
ATOM    592  CA  PRO A 154      36.776  31.077  52.400  1.00 27.08           C  
ATOM    593  C   PRO A 154      37.711  30.586  51.304  1.00 26.26           C  
ATOM    594  O   PRO A 154      38.319  29.526  51.451  1.00 25.82           O  
ATOM    595  CB  PRO A 154      37.600  31.534  53.616  1.00 26.96           C  
ATOM    596  CG  PRO A 154      37.515  30.400  54.587  1.00 27.49           C  
ATOM    597  CD  PRO A 154      36.160  29.797  54.389  1.00 27.46           C  
ATOM    598  N   VAL A 155      37.807  31.363  50.219  1.00 25.90           N  
ATOM    599  CA  VAL A 155      38.855  31.203  49.194  1.00 25.33           C  
ATOM    600  C   VAL A 155      40.202  31.596  49.811  1.00 25.34           C  
ATOM    601  O   VAL A 155      40.234  32.220  50.878  1.00 25.17           O  
ATOM    602  CB  VAL A 155      38.577  32.060  47.913  1.00 25.27           C  
ATOM    603  CG1 VAL A 155      37.234  31.729  47.329  1.00 24.99           C  
ATOM    604  CG2 VAL A 155      38.646  33.566  48.204  1.00 24.78           C  
HETATM  605  N   OCS A 156      41.305  31.256  49.147  1.00 25.48           N  
HETATM  606  CA  OCS A 156      42.632  31.384  49.771  1.00 25.95           C  
HETATM  607  CB  OCS A 156      43.731  30.711  48.952  1.00 25.98           C  
HETATM  608  SG  OCS A 156      43.829  31.171  47.206  1.00 27.12           S  
HETATM  609  C   OCS A 156      43.033  32.804  50.112  1.00 26.25           C  
HETATM  610  O   OCS A 156      43.713  33.028  51.115  1.00 26.44           O  
HETATM  611  OD1 OCS A 156      44.325  32.509  47.083  1.00 24.60           O  
HETATM  612  OD2 OCS A 156      44.703  30.273  46.513  1.00 26.38           O  
HETATM  613  OD3 OCS A 156      42.526  31.031  46.666  1.00 27.05           O  
ATOM    614  N   GLY A 157      42.609  33.758  49.287  1.00 26.33           N  
ATOM    615  CA  GLY A 157      42.918  35.159  49.525  1.00 26.35           C  
ATOM    616  C   GLY A 157      42.189  35.661  50.753  1.00 26.59           C  
ATOM    617  O   GLY A 157      42.801  36.247  51.644  1.00 26.95           O  
ATOM    618  N   LEU A 158      40.881  35.415  50.799  1.00 26.54           N  
ATOM    619  CA  LEU A 158      40.029  35.813  51.917  1.00 26.18           C  
ATOM    620  C   LEU A 158      40.507  35.211  53.247  1.00 26.52           C  
ATOM    621  O   LEU A 158      40.572  35.902  54.263  1.00 26.58           O  
ATOM    622  CB  LEU A 158      38.580  35.411  51.626  1.00 26.07           C  
ATOM    623  CG  LEU A 158      37.475  35.868  52.579  1.00 25.86           C  
ATOM    624  CD1 LEU A 158      37.308  37.370  52.504  1.00 25.15           C  
ATOM    625  CD2 LEU A 158      36.164  35.174  52.235  1.00 25.52           C  
ATOM    626  N   SER A 159      40.859  33.930  53.237  1.00 26.86           N  
ATOM    627  CA  SER A 159      41.368  33.286  54.447  1.00 27.28           C  
ATOM    628  C   SER A 159      42.721  33.865  54.847  1.00 27.17           C  
ATOM    629  O   SER A 159      42.859  34.381  55.948  1.00 27.89           O  
ATOM    630  CB  SER A 159      41.423  31.754  54.304  1.00 27.43           C  
ATOM    631  OG  SER A 159      42.472  31.333  53.438  1.00 28.28           O  
ATOM    632  N   ALA A 160      43.701  33.819  53.949  1.00 27.07           N  
ATOM    633  CA  ALA A 160      45.033  34.340  54.247  1.00 27.01           C  
ATOM    634  C   ALA A 160      44.939  35.756  54.796  1.00 27.21           C  
ATOM    635  O   ALA A 160      45.479  36.036  55.862  1.00 27.42           O  
ATOM    636  CB  ALA A 160      45.940  34.285  53.009  1.00 26.89           C  
ATOM    637  N   ASN A 161      44.229  36.638  54.091  1.00 27.32           N  
ATOM    638  CA  ASN A 161      44.105  38.038  54.516  1.00 27.56           C  
ATOM    639  C   ASN A 161      43.516  38.189  55.911  1.00 27.62           C  
ATOM    640  O   ASN A 161      44.023  38.952  56.722  1.00 27.85           O  
ATOM    641  CB  ASN A 161      43.256  38.862  53.534  1.00 27.70           C  
ATOM    642  CG  ASN A 161      43.931  39.070  52.176  1.00 27.77           C  
ATOM    643  OD1 ASN A 161      45.144  38.941  52.035  1.00 28.84           O  
ATOM    644  ND2 ASN A 161      43.133  39.401  51.175  1.00 27.29           N  
ATOM    645  N   ASN A 162      42.437  37.467  56.185  1.00 27.68           N  
ATOM    646  CA  ASN A 162      41.707  37.658  57.434  1.00 27.43           C  
ATOM    647  C   ASN A 162      42.438  37.091  58.631  1.00 27.20           C  
ATOM    648  O   ASN A 162      42.490  37.714  59.691  1.00 27.31           O  
ATOM    649  CB  ASN A 162      40.293  37.103  57.325  1.00 27.09           C  
ATOM    650  CG  ASN A 162      39.357  38.078  56.651  1.00 28.09           C  
ATOM    651  OD1 ASN A 162      39.138  39.175  57.155  1.00 29.12           O  
ATOM    652  ND2 ASN A 162      38.804  37.692  55.503  1.00 28.44           N  
ATOM    653  N   ILE A 163      43.027  35.919  58.449  1.00 26.76           N  
ATOM    654  CA  ILE A 163      43.792  35.308  59.519  1.00 26.58           C  
ATOM    655  C   ILE A 163      45.071  36.105  59.775  1.00 26.17           C  
ATOM    656  O   ILE A 163      45.509  36.193  60.908  1.00 26.41           O  
ATOM    657  CB  ILE A 163      44.054  33.798  59.262  1.00 26.67           C  
ATOM    658  CG1 ILE A 163      42.730  33.028  59.336  1.00 26.54           C  
ATOM    659  CG2 ILE A 163      45.046  33.229  60.285  1.00 26.83           C  
ATOM    660  CD1 ILE A 163      42.831  31.549  58.982  1.00 26.65           C  
ATOM    661  N   PHE A 164      45.631  36.712  58.729  1.00 25.78           N  
ATOM    662  CA  PHE A 164      46.773  37.621  58.861  1.00 25.45           C  
ATOM    663  C   PHE A 164      46.389  38.810  59.735  1.00 25.77           C  
ATOM    664  O   PHE A 164      47.125  39.169  60.661  1.00 25.57           O  
ATOM    665  CB  PHE A 164      47.245  38.108  57.481  1.00 25.02           C  
ATOM    666  CG  PHE A 164      48.473  38.993  57.522  1.00 24.25           C  
ATOM    667  CD1 PHE A 164      49.744  38.455  57.326  1.00 23.80           C  
ATOM    668  CD2 PHE A 164      48.358  40.363  57.737  1.00 23.87           C  
ATOM    669  CE1 PHE A 164      50.888  39.259  57.357  1.00 23.06           C  
ATOM    670  CE2 PHE A 164      49.494  41.176  57.765  1.00 24.07           C  
ATOM    671  CZ  PHE A 164      50.765  40.615  57.579  1.00 23.98           C  
ATOM    672  N   LYS A 165      45.237  39.410  59.420  1.00 25.96           N  
ATOM    673  CA  LYS A 165      44.670  40.508  60.186  1.00 26.00           C  
ATOM    674  C   LYS A 165      44.529  40.130  61.650  1.00 26.52           C  
ATOM    675  O   LYS A 165      45.015  40.838  62.525  1.00 26.54           O  
ATOM    676  CB  LYS A 165      43.295  40.904  59.630  1.00 26.17           C  
ATOM    677  CG  LYS A 165      43.321  41.976  58.552  1.00 25.55           C  
ATOM    678  CD  LYS A 165      41.914  42.373  58.108  1.00 25.38           C  
ATOM    679  CE  LYS A 165      41.271  43.439  59.012  1.00 24.71           C  
ATOM    680  NZ  LYS A 165      41.860  44.807  58.833  1.00 22.66           N  
ATOM    681  N   LEU A 166      43.866  39.006  61.909  1.00 27.23           N  
ATOM    682  CA  LEU A 166      43.599  38.545  63.276  1.00 27.87           C  
ATOM    683  C   LEU A 166      44.882  38.367  64.114  1.00 28.67           C  
ATOM    684  O   LEU A 166      44.909  38.678  65.310  1.00 28.76           O  
ATOM    685  CB  LEU A 166      42.811  37.235  63.236  1.00 27.53           C  
ATOM    686  CG  LEU A 166      42.231  36.691  64.541  1.00 27.31           C  
ATOM    687  CD1 LEU A 166      41.066  37.522  65.039  1.00 27.71           C  
ATOM    688  CD2 LEU A 166      41.776  35.294  64.317  1.00 27.24           C  
HETATM  689  N   MSE A 167      45.938  37.875  63.477  1.00 28.96           N  
HETATM  690  CA  MSE A 167      47.189  37.605  64.160  1.00 30.25           C  
HETATM  691  C   MSE A 167      48.024  38.847  64.468  1.00 28.80           C  
HETATM  692  O   MSE A 167      48.735  38.886  65.476  1.00 28.70           O  
HETATM  693  CB  MSE A 167      48.014  36.639  63.331  1.00 29.78           C  
HETATM  694  CG  MSE A 167      47.370  35.288  63.191  1.00 31.80           C  
HETATM  695 SE   MSE A 167      48.539  34.058  62.279  1.00 36.57          SE  
HETATM  696  CE  MSE A 167      49.916  33.844  63.654  1.00 34.63           C  
ATOM    697  N   THR A 168      47.927  39.856  63.608  1.00 27.67           N  
ATOM    698  CA  THR A 168      48.832  40.994  63.655  1.00 26.59           C  
ATOM    699  C   THR A 168      48.212  42.222  64.316  1.00 26.19           C  
ATOM    700  O   THR A 168      48.926  43.061  64.852  1.00 25.76           O  
ATOM    701  CB  THR A 168      49.351  41.366  62.237  1.00 26.54           C  
ATOM    702  OG1 THR A 168      48.246  41.499  61.336  1.00 26.55           O  
ATOM    703  CG2 THR A 168      50.276  40.295  61.695  1.00 26.56           C  
ATOM    704  N   GLU A 169      46.886  42.319  64.289  1.00 25.88           N  
ATOM    705  CA  GLU A 169      46.194  43.529  64.717  1.00 25.84           C  
ATOM    706  C   GLU A 169      46.090  43.680  66.218  1.00 25.72           C  
ATOM    707  O   GLU A 169      45.810  42.714  66.921  1.00 26.04           O  
ATOM    708  CB  GLU A 169      44.792  43.578  64.131  1.00 25.54           C  
ATOM    709  CG  GLU A 169      44.770  43.841  62.656  1.00 27.01           C  
ATOM    710  CD  GLU A 169      43.391  44.230  62.131  1.00 28.70           C  
ATOM    711  OE1 GLU A 169      43.347  45.053  61.186  1.00 29.97           O  
ATOM    712  OE2 GLU A 169      42.366  43.728  62.652  1.00 28.48           O  
ATOM    713  N   LYS A 170      46.313  44.899  66.701  1.00 25.68           N  
ATOM    714  CA  LYS A 170      45.905  45.267  68.053  1.00 26.21           C  
ATOM    715  C   LYS A 170      44.823  46.349  68.014  1.00 26.22           C  
ATOM    716  O   LYS A 170      43.636  46.062  68.176  1.00 25.93           O  
ATOM    717  CB  LYS A 170      47.087  45.713  68.930  1.00 26.14           C  
ATOM    718  CG  LYS A 170      46.670  45.911  70.388  1.00 26.24           C  
ATOM    719  CD  LYS A 170      47.801  46.296  71.311  1.00 26.87           C  
ATOM    720  CE  LYS A 170      47.254  46.556  72.714  1.00 27.93           C  
ATOM    721  NZ  LYS A 170      48.127  47.445  73.532  1.00 28.82           N  
ATOM    722  N   ASP A 171      45.253  47.586  67.784  1.00 26.72           N  
ATOM    723  CA  ASP A 171      44.367  48.750  67.797  1.00 27.13           C  
ATOM    724  C   ASP A 171      44.463  49.599  66.521  1.00 26.90           C  
ATOM    725  O   ASP A 171      43.703  50.554  66.351  1.00 26.89           O  
ATOM    726  CB  ASP A 171      44.625  49.604  69.050  1.00 27.28           C  
ATOM    727  CG  ASP A 171      46.109  49.798  69.342  1.00 28.47           C  
ATOM    728  OD1 ASP A 171      46.958  49.440  68.490  1.00 29.56           O  
ATOM    729  OD2 ASP A 171      46.429  50.312  70.436  1.00 29.86           O  
ATOM    730  N   VAL A 172      45.388  49.233  65.633  1.00 26.88           N  
ATOM    731  CA  VAL A 172      45.602  49.925  64.353  1.00 26.67           C  
ATOM    732  C   VAL A 172      45.281  48.958  63.211  1.00 26.68           C  
ATOM    733  O   VAL A 172      45.948  47.925  63.074  1.00 26.44           O  
ATOM    734  CB  VAL A 172      47.078  50.421  64.195  1.00 26.58           C  
ATOM    735  CG1 VAL A 172      47.274  51.173  62.879  1.00 26.20           C  
ATOM    736  CG2 VAL A 172      47.481  51.300  65.353  1.00 26.48           C  
ATOM    737  N   PRO A 173      44.265  49.289  62.386  1.00 26.65           N  
ATOM    738  CA  PRO A 173      43.890  48.431  61.265  1.00 26.68           C  
ATOM    739  C   PRO A 173      45.074  48.140  60.360  1.00 27.00           C  
ATOM    740  O   PRO A 173      45.847  49.045  60.046  1.00 27.15           O  
ATOM    741  CB  PRO A 173      42.859  49.270  60.511  1.00 26.32           C  
ATOM    742  CG  PRO A 173      42.269  50.132  61.534  1.00 26.58           C  
ATOM    743  CD  PRO A 173      43.402  50.481  62.464  1.00 26.64           C  
ATOM    744  N   ILE A 174      45.229  46.873  59.982  1.00 27.31           N  
ATOM    745  CA  ILE A 174      46.216  46.476  58.989  1.00 27.41           C  
ATOM    746  C   ILE A 174      45.483  45.957  57.756  1.00 27.72           C  
ATOM    747  O   ILE A 174      44.655  45.048  57.851  1.00 27.88           O  
ATOM    748  CB  ILE A 174      47.189  45.396  59.535  1.00 27.35           C  
ATOM    749  CG1 ILE A 174      47.968  45.931  60.737  1.00 27.18           C  
ATOM    750  CG2 ILE A 174      48.166  44.940  58.448  1.00 26.80           C  
ATOM    751  CD1 ILE A 174      48.461  44.852  61.685  1.00 26.38           C  
ATOM    752  N   ASP A 175      45.780  46.554  56.608  1.00 28.01           N  
ATOM    753  CA  ASP A 175      45.233  46.107  55.334  1.00 28.40           C  
ATOM    754  C   ASP A 175      46.252  45.222  54.631  1.00 28.50           C  
ATOM    755  O   ASP A 175      47.220  45.722  54.049  1.00 28.50           O  
ATOM    756  CB  ASP A 175      44.864  47.311  54.460  1.00 28.48           C  
ATOM    757  CG  ASP A 175      44.178  46.914  53.166  1.00 29.13           C  
ATOM    758  OD1 ASP A 175      44.074  45.702  52.874  1.00 29.24           O  
ATOM    759  OD2 ASP A 175      43.753  47.830  52.428  1.00 29.88           O  
ATOM    760  N   PRO A 176      46.024  43.898  54.651  1.00 28.81           N  
ATOM    761  CA  PRO A 176      47.020  42.978  54.104  1.00 28.90           C  
ATOM    762  C   PRO A 176      47.161  43.108  52.595  1.00 29.13           C  
ATOM    763  O   PRO A 176      48.158  42.652  52.045  1.00 29.35           O  
ATOM    764  CB  PRO A 176      46.468  41.595  54.464  1.00 28.76           C  
ATOM    765  CG  PRO A 176      45.375  41.832  55.431  1.00 28.81           C  
ATOM    766  CD  PRO A 176      44.831  43.180  55.124  1.00 28.74           C  
ATOM    767  N   THR A 177      46.176  43.722  51.934  1.00 29.66           N  
ATOM    768  CA  THR A 177      46.250  43.958  50.477  1.00 29.98           C  
ATOM    769  C   THR A 177      47.076  45.195  50.143  1.00 30.02           C  
ATOM    770  O   THR A 177      47.527  45.357  49.009  1.00 30.10           O  
ATOM    771  CB  THR A 177      44.853  44.061  49.773  1.00 29.75           C  
ATOM    772  OG1 THR A 177      44.216  45.295  50.116  1.00 30.21           O  
ATOM    773  CG2 THR A 177      43.949  42.895  50.151  1.00 29.75           C  
ATOM    774  N   SER A 178      47.278  46.064  51.128  1.00 30.17           N  
ATOM    775  CA  SER A 178      47.998  47.311  50.884  1.00 30.27           C  
ATOM    776  C   SER A 178      49.413  47.310  51.426  1.00 30.53           C  
ATOM    777  O   SER A 178      50.212  48.162  51.057  1.00 30.68           O  
ATOM    778  CB  SER A 178      47.219  48.501  51.434  1.00 30.15           C  
ATOM    779  OG  SER A 178      45.946  48.585  50.810  1.00 30.30           O  
ATOM    780  N   ILE A 179      49.736  46.363  52.299  1.00 30.95           N  
ATOM    781  CA  ILE A 179      51.092  46.286  52.827  1.00 31.28           C  
ATOM    782  C   ILE A 179      52.099  45.963  51.712  1.00 31.69           C  
ATOM    783  O   ILE A 179      51.759  45.327  50.707  1.00 31.51           O  
ATOM    784  CB  ILE A 179      51.215  45.296  54.011  1.00 31.29           C  
ATOM    785  CG1 ILE A 179      50.850  43.872  53.586  1.00 31.18           C  
ATOM    786  CG2 ILE A 179      50.356  45.764  55.184  1.00 31.46           C  
ATOM    787  CD1 ILE A 179      51.119  42.818  54.644  1.00 31.32           C  
ATOM    788  N   GLU A 180      53.326  46.440  51.894  1.00 32.14           N  
ATOM    789  CA  GLU A 180      54.401  46.245  50.933  1.00 32.66           C  
ATOM    790  C   GLU A 180      54.897  44.794  50.948  1.00 32.72           C  
ATOM    791  O   GLU A 180      55.128  44.225  52.016  1.00 32.66           O  
ATOM    792  CB  GLU A 180      55.547  47.208  51.243  1.00 32.64           C  
ATOM    793  CG  GLU A 180      56.685  47.162  50.246  1.00 34.25           C  
ATOM    794  CD  GLU A 180      57.894  47.954  50.705  1.00 36.71           C  
ATOM    795  OE1 GLU A 180      57.704  49.079  51.227  1.00 37.83           O  
ATOM    796  OE2 GLU A 180      59.035  47.455  50.537  1.00 36.96           O  
ATOM    797  N   TYR A 181      55.053  44.213  49.756  1.00 32.75           N  
ATOM    798  CA  TYR A 181      55.529  42.843  49.596  1.00 32.75           C  
ATOM    799  C   TYR A 181      56.941  42.826  49.044  1.00 33.36           C  
ATOM    800  O   TYR A 181      57.332  43.717  48.293  1.00 33.25           O  
ATOM    801  CB  TYR A 181      54.613  42.057  48.656  1.00 32.18           C  
ATOM    802  CG  TYR A 181      53.299  41.634  49.263  1.00 31.41           C  
ATOM    803  CD1 TYR A 181      53.099  40.327  49.695  1.00 30.57           C  
ATOM    804  CD2 TYR A 181      52.248  42.544  49.403  1.00 31.41           C  
ATOM    805  CE1 TYR A 181      51.879  39.935  50.250  1.00 30.75           C  
ATOM    806  CE2 TYR A 181      51.032  42.168  49.962  1.00 30.57           C  
ATOM    807  CZ  TYR A 181      50.854  40.865  50.380  1.00 30.90           C  
ATOM    808  OH  TYR A 181      49.651  40.501  50.925  1.00 30.58           O  
ATOM    809  N   LEU A 182      57.702  41.804  49.419  1.00 34.39           N  
ATOM    810  CA  LEU A 182      59.062  41.623  48.909  1.00 35.42           C  
ATOM    811  C   LEU A 182      59.047  40.772  47.655  1.00 36.32           C  
ATOM    812  O   LEU A 182      58.323  39.777  47.578  1.00 36.48           O  
ATOM    813  CB  LEU A 182      59.965  40.993  49.969  1.00 35.00           C  
ATOM    814  CG  LEU A 182      60.376  41.932  51.101  1.00 35.20           C  
ATOM    815  CD1 LEU A 182      60.889  41.158  52.315  1.00 34.42           C  
ATOM    816  CD2 LEU A 182      61.414  42.940  50.598  1.00 36.44           C  
ATOM    817  N   GLU A 183      59.855  41.164  46.676  1.00 37.54           N  
ATOM    818  CA  GLU A 183      59.856  40.494  45.381  1.00 38.69           C  
ATOM    819  C   GLU A 183      61.225  39.928  44.972  1.00 39.23           C  
ATOM    820  O   GLU A 183      61.290  38.893  44.298  1.00 39.47           O  
ATOM    821  CB  GLU A 183      59.308  41.438  44.309  1.00 38.64           C  
ATOM    822  CG  GLU A 183      58.816  40.745  43.048  1.00 39.17           C  
ATOM    823  CD  GLU A 183      58.485  41.726  41.941  1.00 39.10           C  
ATOM    824  OE1 GLU A 183      57.521  42.511  42.104  1.00 40.67           O  
ATOM    825  OE2 GLU A 183      59.186  41.705  40.908  1.00 38.92           O  
ATOM    826  N   ASN A 184      62.305  40.599  45.377  1.00 39.90           N  
ATOM    827  CA  ASN A 184      63.669  40.194  44.986  1.00 40.36           C  
ATOM    828  C   ASN A 184      64.518  39.662  46.139  1.00 40.77           C  
ATOM    829  O   ASN A 184      65.748  39.757  46.118  1.00 40.92           O  
ATOM    830  CB  ASN A 184      64.400  41.338  44.267  1.00 40.10           C  
ATOM    831  CG  ASN A 184      64.287  41.248  42.749  1.00 40.27           C  
ATOM    832  OD1 ASN A 184      63.186  41.185  42.193  1.00 40.16           O  
ATOM    833  ND2 ASN A 184      65.434  41.265  42.069  1.00 39.54           N  
ATOM    834  N   THR A 185      63.859  39.115  47.154  1.00 41.27           N  
ATOM    835  CA  THR A 185      64.573  38.447  48.239  1.00 41.90           C  
ATOM    836  C   THR A 185      64.013  37.055  48.488  1.00 41.72           C  
ATOM    837  O   THR A 185      62.799  36.857  48.530  1.00 41.94           O  
ATOM    838  CB  THR A 185      64.594  39.279  49.558  1.00 42.17           C  
ATOM    839  OG1 THR A 185      63.274  39.730  49.880  1.00 42.58           O  
ATOM    840  CG2 THR A 185      65.531  40.498  49.428  1.00 43.11           C  
ATOM    841  N   SER A 186      64.906  36.088  48.643  1.00 41.71           N  
ATOM    842  CA  SER A 186      64.498  34.721  48.915  1.00 41.70           C  
ATOM    843  C   SER A 186      63.984  34.587  50.349  1.00 41.67           C  
ATOM    844  O   SER A 186      64.222  35.456  51.200  1.00 41.69           O  
ATOM    845  CB  SER A 186      65.663  33.759  48.678  1.00 41.71           C  
ATOM    846  OG  SER A 186      66.529  33.725  49.801  1.00 41.75           O  
ATOM    847  N   PHE A 187      63.287  33.489  50.606  1.00 41.43           N  
ATOM    848  CA  PHE A 187      62.749  33.212  51.929  1.00 41.39           C  
ATOM    849  C   PHE A 187      63.843  33.060  52.989  1.00 41.18           C  
ATOM    850  O   PHE A 187      63.660  33.492  54.118  1.00 41.36           O  
ATOM    851  CB  PHE A 187      61.850  31.970  51.876  1.00 41.43           C  
ATOM    852  CG  PHE A 187      61.072  31.720  53.131  1.00 41.29           C  
ATOM    853  CD1 PHE A 187      59.980  32.509  53.458  1.00 41.67           C  
ATOM    854  CD2 PHE A 187      61.423  30.678  53.983  1.00 41.75           C  
ATOM    855  CE1 PHE A 187      59.247  32.268  54.623  1.00 41.91           C  
ATOM    856  CE2 PHE A 187      60.701  30.433  55.150  1.00 41.34           C  
ATOM    857  CZ  PHE A 187      59.611  31.229  55.468  1.00 41.28           C  
ATOM    858  N   ALA A 188      64.972  32.461  52.621  1.00 41.06           N  
ATOM    859  CA  ALA A 188      66.080  32.240  53.562  1.00 41.09           C  
ATOM    860  C   ALA A 188      66.740  33.559  53.950  1.00 41.08           C  
ATOM    861  O   ALA A 188      67.027  33.802  55.121  1.00 41.13           O  
ATOM    862  CB  ALA A 188      67.106  31.272  52.977  1.00 40.87           C  
ATOM    863  N   GLU A 189      66.975  34.404  52.951  1.00 41.15           N  
ATOM    864  CA  GLU A 189      67.446  35.765  53.159  1.00 41.16           C  
ATOM    865  C   GLU A 189      66.669  36.447  54.264  1.00 41.04           C  
ATOM    866  O   GLU A 189      67.256  36.936  55.231  1.00 41.08           O  
ATOM    867  CB  GLU A 189      67.243  36.576  51.890  1.00 41.28           C  
ATOM    868  CG  GLU A 189      68.497  37.028  51.208  1.00 42.47           C  
ATOM    869  CD  GLU A 189      68.210  38.077  50.156  1.00 43.36           C  
ATOM    870  OE1 GLU A 189      68.631  39.239  50.350  1.00 44.14           O  
ATOM    871  OE2 GLU A 189      67.543  37.740  49.151  1.00 43.67           O  
ATOM    872  N   HIS A 190      65.342  36.465  54.111  1.00 40.69           N  
ATOM    873  CA  HIS A 190      64.466  37.229  54.992  1.00 40.12           C  
ATOM    874  C   HIS A 190      64.382  36.672  56.408  1.00 39.95           C  
ATOM    875  O   HIS A 190      64.465  37.424  57.373  1.00 40.10           O  
ATOM    876  CB  HIS A 190      63.070  37.356  54.392  1.00 39.99           C  
ATOM    877  CG  HIS A 190      62.197  38.339  55.109  1.00 40.05           C  
ATOM    878  ND1 HIS A 190      62.493  39.685  55.179  1.00 39.27           N  
ATOM    879  CD2 HIS A 190      61.037  38.171  55.790  1.00 39.42           C  
ATOM    880  CE1 HIS A 190      61.553  40.303  55.873  1.00 39.47           C  
ATOM    881  NE2 HIS A 190      60.657  39.409  56.251  1.00 39.36           N  
ATOM    882  N   VAL A 191      64.217  35.360  56.533  1.00 39.61           N  
ATOM    883  CA  VAL A 191      64.114  34.731  57.845  1.00 39.31           C  
ATOM    884  C   VAL A 191      65.373  35.024  58.670  1.00 39.17           C  
ATOM    885  O   VAL A 191      65.274  35.363  59.851  1.00 39.31           O  
ATOM    886  CB  VAL A 191      63.822  33.209  57.735  1.00 39.36           C  
ATOM    887  CG1 VAL A 191      64.023  32.511  59.068  1.00 39.38           C  
ATOM    888  CG2 VAL A 191      62.400  32.985  57.241  1.00 39.20           C  
ATOM    889  N   ASN A 192      66.539  34.941  58.031  1.00 38.50           N  
ATOM    890  CA  ASN A 192      67.812  35.254  58.677  1.00 38.16           C  
ATOM    891  C   ASN A 192      67.874  36.657  59.298  1.00 37.69           C  
ATOM    892  O   ASN A 192      68.649  36.902  60.227  1.00 37.24           O  
ATOM    893  CB  ASN A 192      68.968  35.075  57.684  1.00 38.37           C  
ATOM    894  CG  ASN A 192      69.328  33.615  57.454  1.00 38.83           C  
ATOM    895  OD1 ASN A 192      68.769  32.710  58.075  1.00 40.17           O  
ATOM    896  ND2 ASN A 192      70.270  33.383  56.555  1.00 39.41           N  
ATOM    897  N   THR A 193      67.047  37.567  58.781  1.00 37.20           N  
ATOM    898  CA  THR A 193      67.061  38.976  59.198  1.00 36.58           C  
ATOM    899  C   THR A 193      66.010  39.285  60.254  1.00 36.31           C  
ATOM    900  O   THR A 193      66.031  40.352  60.868  1.00 36.24           O  
ATOM    901  CB  THR A 193      66.824  39.915  58.013  1.00 36.37           C  
ATOM    902  OG1 THR A 193      65.526  39.665  57.475  1.00 36.05           O  
ATOM    903  CG2 THR A 193      67.882  39.706  56.927  1.00 36.43           C  
ATOM    904  N   LEU A 194      65.089  38.345  60.448  1.00 36.08           N  
ATOM    905  CA  LEU A 194      64.030  38.471  61.437  1.00 35.72           C  
ATOM    906  C   LEU A 194      64.587  38.285  62.830  1.00 35.44           C  
ATOM    907  O   LEU A 194      65.583  37.585  63.005  1.00 35.37           O  
ATOM    908  CB  LEU A 194      62.940  37.442  61.159  1.00 35.83           C  
ATOM    909  CG  LEU A 194      62.230  37.694  59.829  1.00 35.79           C  
ATOM    910  CD1 LEU A 194      61.172  36.653  59.581  1.00 35.86           C  
ATOM    911  CD2 LEU A 194      61.626  39.094  59.803  1.00 35.71           C  
ATOM    912  N   ASP A 195      63.888  38.772  63.818  1.00 35.29           N  
ATOM    913  CA  ASP A 195      64.384  38.762  65.156  1.00 35.37           C  
ATOM    914  C   ASP A 195      63.954  37.526  65.864  1.00 35.14           C  
ATOM    915  O   ASP A 195      62.916  37.009  65.616  1.00 36.18           O  
ATOM    916  CB  ASP A 195      63.833  39.974  65.861  1.00 35.69           C  
ATOM    917  CG  ASP A 195      64.206  40.039  67.284  1.00 35.86           C  
ATOM    918  OD1 ASP A 195      64.343  38.997  67.921  1.00 36.44           O  
ATOM    919  OD2 ASP A 195      64.339  41.152  67.768  1.00 36.20           O  
ATOM    920  N   SER A 196      64.764  37.051  66.771  1.00 34.44           N  
ATOM    921  CA  SER A 196      64.646  35.689  67.241  1.00 33.75           C  
ATOM    922  C   SER A 196      63.600  35.500  68.311  1.00 32.85           C  
ATOM    923  O   SER A 196      62.946  34.504  68.361  1.00 32.76           O  
ATOM    924  CB  SER A 196      65.996  35.206  67.711  1.00 33.88           C  
ATOM    925  OG  SER A 196      66.986  35.884  66.961  1.00 34.58           O  
ATOM    926  N   HIS A 197      63.457  36.476  69.166  1.00 31.67           N  
ATOM    927  CA  HIS A 197      62.525  36.407  70.260  1.00 30.84           C  
ATOM    928  C   HIS A 197      61.216  37.090  69.875  1.00 29.75           C  
ATOM    929  O   HIS A 197      60.372  37.347  70.727  1.00 29.61           O  
ATOM    930  CB  HIS A 197      63.155  36.997  71.529  1.00 31.23           C  
ATOM    931  CG  HIS A 197      64.477  36.380  71.893  1.00 32.25           C  
ATOM    932  ND1 HIS A 197      65.636  37.121  72.029  1.00 33.35           N  
ATOM    933  CD2 HIS A 197      64.825  35.090  72.133  1.00 33.01           C  
ATOM    934  CE1 HIS A 197      66.637  36.316  72.345  1.00 33.42           C  
ATOM    935  NE2 HIS A 197      66.172  35.079  72.413  1.00 33.42           N  
ATOM    936  N   LYS A 198      61.042  37.349  68.583  1.00 28.44           N  
ATOM    937  CA  LYS A 198      59.792  37.903  68.072  1.00 27.69           C  
ATOM    938  C   LYS A 198      58.973  36.885  67.262  1.00 26.83           C  
ATOM    939  O   LYS A 198      59.531  35.961  66.683  1.00 26.66           O  
ATOM    940  CB  LYS A 198      60.071  39.155  67.239  1.00 27.93           C  
ATOM    941  CG  LYS A 198      60.687  40.295  68.036  1.00 28.62           C  
ATOM    942  CD  LYS A 198      60.142  41.655  67.606  1.00 30.48           C  
ATOM    943  CE  LYS A 198      60.744  42.799  68.429  1.00 31.46           C  
ATOM    944  NZ  LYS A 198      60.709  42.557  69.916  1.00 31.38           N  
ATOM    945  N   ASN A 199      57.650  37.050  67.239  1.00 25.98           N  
ATOM    946  CA  ASN A 199      56.777  36.222  66.395  1.00 25.55           C  
ATOM    947  C   ASN A 199      56.496  36.946  65.095  1.00 25.15           C  
ATOM    948  O   ASN A 199      56.274  38.154  65.103  1.00 25.37           O  
ATOM    949  CB  ASN A 199      55.420  35.926  67.055  1.00 25.59           C  
ATOM    950  CG  ASN A 199      55.530  35.530  68.514  1.00 25.71           C  
ATOM    951  OD1 ASN A 199      56.067  34.472  68.855  1.00 25.95           O  
ATOM    952  ND2 ASN A 199      54.985  36.372  69.388  1.00 26.13           N  
ATOM    953  N   TYR A 200      56.496  36.211  63.986  1.00 24.53           N  
ATOM    954  CA  TYR A 200      56.150  36.771  62.686  1.00 24.03           C  
ATOM    955  C   TYR A 200      55.170  35.869  61.991  1.00 24.05           C  
ATOM    956  O   TYR A 200      55.192  34.650  62.181  1.00 23.99           O  
ATOM    957  CB  TYR A 200      57.374  36.907  61.787  1.00 23.66           C  
ATOM    958  CG  TYR A 200      58.449  37.785  62.348  1.00 23.80           C  
ATOM    959  CD1 TYR A 200      58.462  39.156  62.085  1.00 23.56           C  
ATOM    960  CD2 TYR A 200      59.463  37.251  63.144  1.00 22.97           C  
ATOM    961  CE1 TYR A 200      59.452  39.965  62.603  1.00 23.04           C  
ATOM    962  CE2 TYR A 200      60.454  38.059  63.665  1.00 23.20           C  
ATOM    963  CZ  TYR A 200      60.444  39.407  63.385  1.00 22.62           C  
ATOM    964  OH  TYR A 200      61.417  40.204  63.902  1.00 23.33           O  
ATOM    965  N   VAL A 201      54.308  36.473  61.181  1.00 23.93           N  
ATOM    966  CA  VAL A 201      53.517  35.721  60.227  1.00 23.87           C  
ATOM    967  C   VAL A 201      53.844  36.257  58.833  1.00 23.74           C  
ATOM    968  O   VAL A 201      53.940  37.466  58.638  1.00 23.74           O  
ATOM    969  CB  VAL A 201      52.009  35.774  60.572  1.00 24.04           C  
ATOM    970  CG1 VAL A 201      51.468  37.208  60.482  1.00 23.74           C  
ATOM    971  CG2 VAL A 201      51.215  34.812  59.683  1.00 24.44           C  
ATOM    972  N   VAL A 202      54.064  35.362  57.878  1.00 23.75           N  
ATOM    973  CA  VAL A 202      54.396  35.771  56.516  1.00 23.88           C  
ATOM    974  C   VAL A 202      53.205  35.494  55.616  1.00 24.25           C  
ATOM    975  O   VAL A 202      52.682  34.374  55.610  1.00 24.36           O  
ATOM    976  CB  VAL A 202      55.614  35.000  55.960  1.00 23.92           C  
ATOM    977  CG1 VAL A 202      56.025  35.575  54.621  1.00 23.94           C  
ATOM    978  CG2 VAL A 202      56.786  35.035  56.934  1.00 23.78           C  
ATOM    979  N   ILE A 203      52.762  36.512  54.878  1.00 24.36           N  
ATOM    980  CA  ILE A 203      51.699  36.337  53.882  1.00 24.58           C  
ATOM    981  C   ILE A 203      52.321  36.158  52.505  1.00 24.96           C  
ATOM    982  O   ILE A 203      53.091  37.013  52.049  1.00 25.03           O  
ATOM    983  CB  ILE A 203      50.656  37.495  53.883  1.00 24.39           C  
ATOM    984  CG1 ILE A 203      49.602  37.273  52.786  1.00 23.73           C  
ATOM    985  CG2 ILE A 203      51.338  38.858  53.741  1.00 24.66           C  
ATOM    986  CD1 ILE A 203      48.185  37.737  53.157  1.00 22.83           C  
ATOM    987  N   VAL A 204      51.987  35.040  51.860  1.00 25.44           N  
ATOM    988  CA  VAL A 204      52.602  34.649  50.588  1.00 25.64           C  
ATOM    989  C   VAL A 204      51.612  34.687  49.429  1.00 26.04           C  
ATOM    990  O   VAL A 204      50.552  34.046  49.466  1.00 25.98           O  
ATOM    991  CB  VAL A 204      53.259  33.255  50.663  1.00 25.70           C  
ATOM    992  CG1 VAL A 204      54.038  32.946  49.379  1.00 25.11           C  
ATOM    993  CG2 VAL A 204      54.183  33.170  51.876  1.00 25.98           C  
ATOM    994  N   ASN A 205      51.977  35.456  48.407  1.00 26.16           N  
ATOM    995  CA  ASN A 205      51.280  35.448  47.146  1.00 26.51           C  
ATOM    996  C   ASN A 205      52.088  34.682  46.096  1.00 26.79           C  
ATOM    997  O   ASN A 205      53.024  35.225  45.496  1.00 26.81           O  
ATOM    998  CB  ASN A 205      51.021  36.882  46.695  1.00 26.55           C  
ATOM    999  CG  ASN A 205      50.028  36.967  45.560  1.00 26.51           C  
ATOM   1000  OD1 ASN A 205      49.920  36.065  44.735  1.00 26.64           O  
ATOM   1001  ND2 ASN A 205      49.295  38.061  45.511  1.00 27.96           N  
ATOM   1002  N   ASP A 206      51.720  33.424  45.869  1.00 27.15           N  
ATOM   1003  CA  ASP A 206      52.453  32.580  44.924  1.00 27.90           C  
ATOM   1004  C   ASP A 206      51.851  32.552  43.509  1.00 28.20           C  
ATOM   1005  O   ASP A 206      50.901  31.810  43.238  1.00 28.04           O  
ATOM   1006  CB  ASP A 206      52.636  31.157  45.472  1.00 27.92           C  
ATOM   1007  CG  ASP A 206      53.772  30.410  44.784  1.00 28.42           C  
ATOM   1008  OD1 ASP A 206      53.807  30.397  43.532  1.00 28.86           O  
ATOM   1009  OD2 ASP A 206      54.636  29.839  45.493  1.00 28.07           O  
ATOM   1010  N   GLY A 207      52.441  33.343  42.612  1.00 28.54           N  
ATOM   1011  CA  GLY A 207      52.012  33.414  41.219  1.00 29.16           C  
ATOM   1012  C   GLY A 207      52.165  32.122  40.432  1.00 29.77           C  
ATOM   1013  O   GLY A 207      51.548  31.956  39.377  1.00 30.23           O  
ATOM   1014  N   ARG A 208      52.983  31.203  40.929  1.00 29.90           N  
ATOM   1015  CA  ARG A 208      53.179  29.928  40.243  1.00 30.17           C  
ATOM   1016  C   ARG A 208      52.020  28.982  40.534  1.00 30.22           C  
ATOM   1017  O   ARG A 208      51.553  28.265  39.646  1.00 30.46           O  
ATOM   1018  CB  ARG A 208      54.523  29.292  40.625  1.00 30.17           C  
ATOM   1019  CG  ARG A 208      55.736  30.075  40.138  1.00 30.65           C  
ATOM   1020  CD  ARG A 208      56.167  31.133  41.137  1.00 31.32           C  
ATOM   1021  NE  ARG A 208      57.409  30.752  41.796  1.00 32.61           N  
ATOM   1022  CZ  ARG A 208      57.493  30.084  42.942  1.00 32.54           C  
ATOM   1023  NH1 ARG A 208      56.399  29.714  43.591  1.00 32.58           N  
ATOM   1024  NH2 ARG A 208      58.684  29.788  43.439  1.00 31.90           N  
ATOM   1025  N   LEU A 209      51.563  28.997  41.784  1.00 29.99           N  
ATOM   1026  CA  LEU A 209      50.438  28.188  42.224  1.00 29.55           C  
ATOM   1027  C   LEU A 209      49.110  28.936  42.048  1.00 29.45           C  
ATOM   1028  O   LEU A 209      48.041  28.324  42.077  1.00 29.63           O  
ATOM   1029  CB  LEU A 209      50.619  27.833  43.698  1.00 29.66           C  
ATOM   1030  CG  LEU A 209      51.639  26.799  44.195  1.00 30.02           C  
ATOM   1031  CD1 LEU A 209      53.078  27.071  43.753  1.00 31.40           C  
ATOM   1032  CD2 LEU A 209      51.582  26.746  45.706  1.00 29.55           C  
ATOM   1033  N   GLY A 210      49.181  30.256  41.865  1.00 29.03           N  
ATOM   1034  CA  GLY A 210      47.996  31.116  41.898  1.00 28.55           C  
ATOM   1035  C   GLY A 210      47.318  31.017  43.255  1.00 28.21           C  
ATOM   1036  O   GLY A 210      46.103  30.808  43.344  1.00 28.24           O  
ATOM   1037  N   HIS A 211      48.109  31.172  44.313  1.00 27.65           N  
ATOM   1038  CA  HIS A 211      47.678  30.763  45.641  1.00 27.29           C  
ATOM   1039  C   HIS A 211      48.259  31.637  46.720  1.00 26.97           C  
ATOM   1040  O   HIS A 211      49.465  31.871  46.756  1.00 26.86           O  
ATOM   1041  CB  HIS A 211      48.070  29.296  45.892  1.00 27.11           C  
ATOM   1042  CG  HIS A 211      47.586  28.738  47.200  1.00 26.90           C  
ATOM   1043  ND1 HIS A 211      46.267  28.803  47.604  1.00 25.76           N  
ATOM   1044  CD2 HIS A 211      48.243  28.069  48.177  1.00 25.18           C  
ATOM   1045  CE1 HIS A 211      46.137  28.213  48.778  1.00 24.08           C  
ATOM   1046  NE2 HIS A 211      47.321  27.761  49.149  1.00 24.66           N  
ATOM   1047  N   LYS A 212      47.376  32.116  47.593  1.00 26.99           N  
ATOM   1048  CA  LYS A 212      47.763  32.807  48.816  1.00 26.58           C  
ATOM   1049  C   LYS A 212      47.687  31.835  49.984  1.00 25.94           C  
ATOM   1050  O   LYS A 212      46.807  30.977  50.024  1.00 25.88           O  
ATOM   1051  CB  LYS A 212      46.851  34.002  49.079  1.00 27.01           C  
ATOM   1052  CG  LYS A 212      47.215  35.260  48.292  1.00 27.88           C  
ATOM   1053  CD  LYS A 212      46.663  36.510  48.984  1.00 29.69           C  
ATOM   1054  CE  LYS A 212      46.915  37.797  48.178  1.00 30.28           C  
ATOM   1055  NZ  LYS A 212      46.172  37.820  46.874  1.00 31.19           N  
ATOM   1056  N   PHE A 213      48.622  31.976  50.919  1.00 25.10           N  
ATOM   1057  CA  PHE A 213      48.655  31.203  52.158  1.00 24.33           C  
ATOM   1058  C   PHE A 213      49.554  31.910  53.173  1.00 24.04           C  
ATOM   1059  O   PHE A 213      50.342  32.797  52.820  1.00 23.59           O  
ATOM   1060  CB  PHE A 213      49.126  29.749  51.918  1.00 24.31           C  
ATOM   1061  CG  PHE A 213      50.533  29.629  51.375  1.00 24.42           C  
ATOM   1062  CD1 PHE A 213      51.609  29.377  52.236  1.00 25.16           C  
ATOM   1063  CD2 PHE A 213      50.784  29.751  50.003  1.00 24.04           C  
ATOM   1064  CE1 PHE A 213      52.928  29.249  51.746  1.00 24.65           C  
ATOM   1065  CE2 PHE A 213      52.081  29.627  49.491  1.00 24.77           C  
ATOM   1066  CZ  PHE A 213      53.166  29.370  50.367  1.00 25.33           C  
ATOM   1067  N   LEU A 214      49.431  31.514  54.437  1.00 23.64           N  
ATOM   1068  CA  LEU A 214      50.248  32.091  55.501  1.00 22.76           C  
ATOM   1069  C   LEU A 214      51.320  31.124  55.992  1.00 22.15           C  
ATOM   1070  O   LEU A 214      51.121  29.907  56.002  1.00 22.23           O  
ATOM   1071  CB  LEU A 214      49.368  32.529  56.666  1.00 22.49           C  
ATOM   1072  CG  LEU A 214      48.218  33.496  56.404  1.00 22.70           C  
ATOM   1073  CD1 LEU A 214      47.480  33.752  57.699  1.00 22.01           C  
ATOM   1074  CD2 LEU A 214      48.686  34.820  55.790  1.00 23.45           C  
ATOM   1075  N   ILE A 215      52.467  31.674  56.373  1.00 21.65           N  
ATOM   1076  CA  ILE A 215      53.500  30.914  57.054  1.00 21.14           C  
ATOM   1077  C   ILE A 215      53.665  31.537  58.430  1.00 21.37           C  
ATOM   1078  O   ILE A 215      54.056  32.698  58.547  1.00 21.60           O  
ATOM   1079  CB  ILE A 215      54.850  30.944  56.312  1.00 21.00           C  
ATOM   1080  CG1 ILE A 215      54.686  30.437  54.879  1.00 20.22           C  
ATOM   1081  CG2 ILE A 215      55.905  30.131  57.087  1.00 20.34           C  
ATOM   1082  CD1 ILE A 215      55.914  30.630  54.014  1.00 20.76           C  
ATOM   1083  N   ASP A 216      53.355  30.765  59.467  1.00 21.45           N  
ATOM   1084  CA  ASP A 216      53.449  31.239  60.832  1.00 21.79           C  
ATOM   1085  C   ASP A 216      54.826  30.918  61.421  1.00 22.00           C  
ATOM   1086  O   ASP A 216      55.294  29.781  61.357  1.00 22.20           O  
ATOM   1087  CB  ASP A 216      52.334  30.620  61.668  1.00 21.80           C  
ATOM   1088  CG  ASP A 216      52.248  31.208  63.062  1.00 23.11           C  
ATOM   1089  OD1 ASP A 216      53.067  32.088  63.415  1.00 24.35           O  
ATOM   1090  OD2 ASP A 216      51.344  30.789  63.813  1.00 25.23           O  
ATOM   1091  N   LEU A 217      55.475  31.925  61.990  1.00 21.95           N  
ATOM   1092  CA  LEU A 217      56.802  31.734  62.564  1.00 22.19           C  
ATOM   1093  C   LEU A 217      56.859  32.215  64.009  1.00 22.28           C  
ATOM   1094  O   LEU A 217      57.446  33.259  64.284  1.00 22.68           O  
ATOM   1095  CB  LEU A 217      57.855  32.476  61.736  1.00 22.00           C  
ATOM   1096  CG  LEU A 217      57.998  32.138  60.251  1.00 22.13           C  
ATOM   1097  CD1 LEU A 217      58.840  33.214  59.596  1.00 23.00           C  
ATOM   1098  CD2 LEU A 217      58.629  30.763  60.038  1.00 22.04           C  
ATOM   1099  N   PRO A 218      56.243  31.466  64.942  1.00 22.20           N  
ATOM   1100  CA  PRO A 218      56.339  31.872  66.342  1.00 22.16           C  
ATOM   1101  C   PRO A 218      57.751  31.680  66.863  1.00 22.31           C  
ATOM   1102  O   PRO A 218      58.484  30.849  66.334  1.00 22.93           O  
ATOM   1103  CB  PRO A 218      55.387  30.912  67.047  1.00 21.95           C  
ATOM   1104  CG  PRO A 218      55.299  29.737  66.166  1.00 21.80           C  
ATOM   1105  CD  PRO A 218      55.443  30.240  64.776  1.00 22.04           C  
ATOM   1106  N   ALA A 219      58.137  32.448  67.874  1.00 22.43           N  
ATOM   1107  CA  ALA A 219      59.447  32.268  68.519  1.00 22.68           C  
ATOM   1108  C   ALA A 219      59.474  31.034  69.437  1.00 22.55           C  
ATOM   1109  O   ALA A 219      58.609  30.869  70.308  1.00 22.44           O  
ATOM   1110  CB  ALA A 219      59.839  33.520  69.292  1.00 22.49           C  
ATOM   1111  N   ARG A 225      64.582  27.240  64.207  1.00 30.49           N  
ATOM   1112  CA  ARG A 225      63.283  27.896  64.071  1.00 30.59           C  
ATOM   1113  C   ARG A 225      62.245  26.979  63.422  1.00 30.30           C  
ATOM   1114  O   ARG A 225      62.591  26.054  62.678  1.00 30.65           O  
ATOM   1115  CB  ARG A 225      63.405  29.190  63.267  1.00 30.53           C  
ATOM   1116  CG  ARG A 225      62.663  30.317  63.925  1.00 31.45           C  
ATOM   1117  CD  ARG A 225      62.296  31.452  62.993  1.00 32.04           C  
ATOM   1118  NE  ARG A 225      61.163  32.165  63.586  1.00 32.26           N  
ATOM   1119  CZ  ARG A 225      61.232  33.356  64.169  1.00 31.00           C  
ATOM   1120  NH1 ARG A 225      62.377  34.024  64.215  1.00 30.71           N  
ATOM   1121  NH2 ARG A 225      60.136  33.889  64.690  1.00 30.74           N  
ATOM   1122  N   THR A 226      60.973  27.231  63.709  1.00 29.43           N  
ATOM   1123  CA  THR A 226      59.923  26.345  63.238  1.00 28.90           C  
ATOM   1124  C   THR A 226      58.778  27.115  62.591  1.00 27.92           C  
ATOM   1125  O   THR A 226      58.381  28.167  63.078  1.00 27.80           O  
ATOM   1126  CB  THR A 226      59.372  25.443  64.372  1.00 29.24           C  
ATOM   1127  OG1 THR A 226      58.403  26.164  65.132  1.00 30.61           O  
ATOM   1128  CG2 THR A 226      60.487  24.952  65.308  1.00 29.51           C  
ATOM   1129  N   ALA A 227      58.251  26.560  61.499  1.00 27.06           N  
ATOM   1130  CA  ALA A 227      57.222  27.198  60.678  1.00 25.79           C  
ATOM   1131  C   ALA A 227      55.966  26.340  60.563  1.00 25.25           C  
ATOM   1132  O   ALA A 227      56.039  25.114  60.657  1.00 25.23           O  
ATOM   1133  CB  ALA A 227      57.775  27.460  59.307  1.00 25.79           C  
ATOM   1134  N   TYR A 228      54.819  26.984  60.346  1.00 24.71           N  
ATOM   1135  CA  TYR A 228      53.542  26.287  60.126  1.00 24.25           C  
ATOM   1136  C   TYR A 228      52.774  26.940  59.004  1.00 24.14           C  
ATOM   1137  O   TYR A 228      52.905  28.130  58.763  1.00 24.04           O  
ATOM   1138  CB  TYR A 228      52.647  26.304  61.367  1.00 24.16           C  
ATOM   1139  CG  TYR A 228      53.288  25.800  62.627  1.00 24.13           C  
ATOM   1140  CD1 TYR A 228      53.155  24.472  63.014  1.00 24.59           C  
ATOM   1141  CD2 TYR A 228      54.016  26.659  63.447  1.00 24.83           C  
ATOM   1142  CE1 TYR A 228      53.749  24.001  64.191  1.00 24.91           C  
ATOM   1143  CE2 TYR A 228      54.616  26.200  64.617  1.00 25.52           C  
ATOM   1144  CZ  TYR A 228      54.472  24.874  64.985  1.00 24.68           C  
ATOM   1145  OH  TYR A 228      55.058  24.433  66.139  1.00 23.93           O  
ATOM   1146  N   ILE A 229      51.942  26.152  58.339  1.00 24.17           N  
ATOM   1147  CA  ILE A 229      51.145  26.641  57.230  1.00 24.07           C  
ATOM   1148  C   ILE A 229      49.716  26.870  57.689  1.00 24.25           C  
ATOM   1149  O   ILE A 229      49.101  25.986  58.291  1.00 24.47           O  
ATOM   1150  CB  ILE A 229      51.145  25.632  56.055  1.00 23.92           C  
ATOM   1151  CG1 ILE A 229      52.562  25.434  55.506  1.00 23.33           C  
ATOM   1152  CG2 ILE A 229      50.122  26.026  54.962  1.00 23.74           C  
ATOM   1153  CD1 ILE A 229      53.229  26.691  55.008  1.00 23.05           C  
ATOM   1154  N   ILE A 230      49.204  28.068  57.427  1.00 24.00           N  
ATOM   1155  CA  ILE A 230      47.777  28.294  57.478  1.00 24.00           C  
ATOM   1156  C   ILE A 230      47.328  28.631  56.062  1.00 24.07           C  
ATOM   1157  O   ILE A 230      47.916  29.509  55.413  1.00 23.84           O  
ATOM   1158  CB  ILE A 230      47.383  29.427  58.451  1.00 24.21           C  
ATOM   1159  CG1 ILE A 230      47.868  29.122  59.865  1.00 23.93           C  
ATOM   1160  CG2 ILE A 230      45.860  29.618  58.474  1.00 24.11           C  
ATOM   1161  CD1 ILE A 230      48.005  30.349  60.735  1.00 24.63           C  
ATOM   1162  N   GLN A 231      46.292  27.928  55.596  1.00 23.95           N  
ATOM   1163  CA  GLN A 231      45.744  28.120  54.257  1.00 23.65           C  
ATOM   1164  C   GLN A 231      44.318  27.601  54.117  1.00 23.74           C  
ATOM   1165  O   GLN A 231      43.844  26.817  54.948  1.00 23.53           O  
ATOM   1166  CB  GLN A 231      46.618  27.422  53.225  1.00 23.81           C  
ATOM   1167  CG  GLN A 231      46.539  25.907  53.253  1.00 24.49           C  
ATOM   1168  CD  GLN A 231      47.158  25.266  52.028  1.00 26.31           C  
ATOM   1169  OE1 GLN A 231      47.405  25.930  51.019  1.00 27.72           O  
ATOM   1170  NE2 GLN A 231      47.413  23.964  52.108  1.00 26.39           N  
ATOM   1171  N   SER A 232      43.642  28.072  53.069  1.00 23.57           N  
ATOM   1172  CA  SER A 232      42.456  27.412  52.521  1.00 23.49           C  
ATOM   1173  C   SER A 232      42.641  27.296  50.998  1.00 23.61           C  
ATOM   1174  O   SER A 232      43.547  27.919  50.427  1.00 23.35           O  
ATOM   1175  CB  SER A 232      41.198  28.212  52.836  1.00 23.61           C  
ATOM   1176  OG  SER A 232      41.117  29.367  52.015  1.00 23.23           O  
ATOM   1177  N   ASP A 233      41.793  26.506  50.340  1.00 23.39           N  
ATOM   1178  CA  ASP A 233      41.882  26.342  48.884  1.00 23.27           C  
ATOM   1179  C   ASP A 233      40.602  25.793  48.261  1.00 22.75           C  
ATOM   1180  O   ASP A 233      40.120  24.737  48.652  1.00 22.81           O  
ATOM   1181  CB  ASP A 233      43.047  25.418  48.524  1.00 23.38           C  
ATOM   1182  CG  ASP A 233      43.626  25.723  47.175  1.00 23.95           C  
ATOM   1183  OD1 ASP A 233      42.863  26.102  46.269  1.00 25.92           O  
ATOM   1184  OD2 ASP A 233      44.853  25.590  47.013  1.00 26.46           O  
ATOM   1185  N   LEU A 234      40.066  26.501  47.279  1.00 21.91           N  
ATOM   1186  CA  LEU A 234      38.909  26.003  46.552  1.00 21.61           C  
ATOM   1187  C   LEU A 234      39.294  25.007  45.451  1.00 21.54           C  
ATOM   1188  O   LEU A 234      38.431  24.390  44.834  1.00 21.74           O  
ATOM   1189  CB  LEU A 234      38.082  27.158  45.987  1.00 21.14           C  
ATOM   1190  CG  LEU A 234      37.303  27.973  47.019  1.00 21.22           C  
ATOM   1191  CD1 LEU A 234      36.099  28.637  46.346  1.00 21.88           C  
ATOM   1192  CD2 LEU A 234      36.839  27.123  48.210  1.00 20.21           C  
ATOM   1193  N   GLY A 235      40.592  24.848  45.225  1.00 21.40           N  
ATOM   1194  CA  GLY A 235      41.097  23.961  44.191  1.00 21.43           C  
ATOM   1195  C   GLY A 235      41.055  24.642  42.841  1.00 21.61           C  
ATOM   1196  O   GLY A 235      40.786  25.839  42.753  1.00 21.80           O  
ATOM   1197  N   GLY A 236      41.311  23.875  41.787  1.00 21.65           N  
ATOM   1198  CA  GLY A 236      41.306  24.409  40.431  1.00 21.41           C  
ATOM   1199  C   GLY A 236      42.707  24.588  39.889  1.00 21.29           C  
ATOM   1200  O   GLY A 236      42.909  24.591  38.674  1.00 21.54           O  
ATOM   1201  N   GLY A 237      43.673  24.746  40.790  1.00 20.84           N  
ATOM   1202  CA  GLY A 237      45.075  24.839  40.405  1.00 20.60           C  
ATOM   1203  C   GLY A 237      45.827  23.552  40.695  1.00 20.47           C  
ATOM   1204  O   GLY A 237      45.239  22.465  40.722  1.00 20.24           O  
ATOM   1205  N   ALA A 238      47.133  23.678  40.923  1.00 20.30           N  
ATOM   1206  CA  ALA A 238      47.992  22.523  41.150  1.00 20.13           C  
ATOM   1207  C   ALA A 238      47.582  21.726  42.393  1.00 20.37           C  
ATOM   1208  O   ALA A 238      47.771  20.509  42.451  1.00 20.23           O  
ATOM   1209  CB  ALA A 238      49.439  22.957  41.244  1.00 19.76           C  
ATOM   1210  N   LEU A 239      47.015  22.414  43.377  1.00 20.61           N  
ATOM   1211  CA  LEU A 239      46.658  21.781  44.636  1.00 21.22           C  
ATOM   1212  C   LEU A 239      45.164  21.491  44.710  1.00 21.30           C  
ATOM   1213  O   LEU A 239      44.356  22.242  44.152  1.00 21.23           O  
ATOM   1214  CB  LEU A 239      47.069  22.659  45.821  1.00 21.39           C  
ATOM   1215  CG  LEU A 239      48.552  22.952  46.026  1.00 22.17           C  
ATOM   1216  CD1 LEU A 239      48.679  24.292  46.721  1.00 23.04           C  
ATOM   1217  CD2 LEU A 239      49.244  21.843  46.818  1.00 22.71           C  
ATOM   1218  N   PRO A 240      44.794  20.405  45.412  1.00 21.29           N  
ATOM   1219  CA  PRO A 240      43.388  20.088  45.635  1.00 21.52           C  
ATOM   1220  C   PRO A 240      42.754  21.049  46.639  1.00 21.53           C  
ATOM   1221  O   PRO A 240      43.466  21.641  47.463  1.00 21.66           O  
ATOM   1222  CB  PRO A 240      43.441  18.689  46.253  1.00 21.46           C  
ATOM   1223  CG  PRO A 240      44.760  18.644  46.950  1.00 21.55           C  
ATOM   1224  CD  PRO A 240      45.687  19.424  46.055  1.00 21.29           C  
ATOM   1225  N   ALA A 241      41.429  21.170  46.572  1.00 21.25           N  
ATOM   1226  CA  ALA A 241      40.645  21.925  47.537  1.00 20.99           C  
ATOM   1227  C   ALA A 241      40.922  21.468  48.965  1.00 21.22           C  
ATOM   1228  O   ALA A 241      41.162  20.284  49.216  1.00 21.48           O  
ATOM   1229  CB  ALA A 241      39.159  21.803  47.218  1.00 20.89           C  
ATOM   1230  N   VAL A 242      40.923  22.421  49.894  1.00 21.38           N  
ATOM   1231  CA  VAL A 242      41.043  22.132  51.331  1.00 21.03           C  
ATOM   1232  C   VAL A 242      40.404  23.260  52.132  1.00 21.15           C  
ATOM   1233  O   VAL A 242      40.757  24.425  51.952  1.00 20.98           O  
ATOM   1234  CB  VAL A 242      42.526  21.845  51.776  1.00 21.03           C  
ATOM   1235  CG1 VAL A 242      43.435  23.083  51.667  1.00 20.14           C  
ATOM   1236  CG2 VAL A 242      42.565  21.247  53.171  1.00 20.50           C  
ATOM   1237  N   ARG A 243      39.426  22.913  52.969  1.00 21.42           N  
ATOM   1238  CA  ARG A 243      38.835  23.875  53.892  1.00 21.76           C  
ATOM   1239  C   ARG A 243      39.857  24.223  54.967  1.00 21.93           C  
ATOM   1240  O   ARG A 243      40.585  23.359  55.463  1.00 21.45           O  
ATOM   1241  CB  ARG A 243      37.570  23.320  54.552  1.00 21.87           C  
ATOM   1242  CG  ARG A 243      36.346  23.268  53.655  1.00 22.18           C  
ATOM   1243  CD  ARG A 243      35.302  22.273  54.215  1.00 23.67           C  
ATOM   1244  NE  ARG A 243      34.586  22.802  55.377  1.00 24.60           N  
ATOM   1245  CZ  ARG A 243      34.801  22.439  56.641  1.00 26.31           C  
ATOM   1246  NH1 ARG A 243      35.725  21.526  56.946  1.00 26.05           N  
ATOM   1247  NH2 ARG A 243      34.088  22.999  57.614  1.00 26.70           N  
ATOM   1248  N   VAL A 244      39.909  25.499  55.318  1.00 22.56           N  
ATOM   1249  CA  VAL A 244      40.789  25.955  56.384  1.00 23.26           C  
ATOM   1250  C   VAL A 244      40.576  25.181  57.702  1.00 23.87           C  
ATOM   1251  O   VAL A 244      41.550  24.863  58.382  1.00 24.61           O  
ATOM   1252  CB  VAL A 244      40.681  27.492  56.589  1.00 23.39           C  
ATOM   1253  CG1 VAL A 244      39.274  27.902  57.114  1.00 22.42           C  
ATOM   1254  CG2 VAL A 244      41.818  28.000  57.484  1.00 23.14           C  
ATOM   1255  N   GLU A 245      39.323  24.858  58.039  1.00 24.12           N  
ATOM   1256  CA  GLU A 245      39.012  24.051  59.231  1.00 24.64           C  
ATOM   1257  C   GLU A 245      39.750  22.718  59.184  1.00 24.96           C  
ATOM   1258  O   GLU A 245      40.441  22.362  60.126  1.00 25.24           O  
ATOM   1259  CB  GLU A 245      37.504  23.797  59.381  1.00 24.30           C  
ATOM   1260  CG  GLU A 245      36.652  25.037  59.625  1.00 25.22           C  
ATOM   1261  CD  GLU A 245      36.272  25.781  58.350  1.00 26.43           C  
ATOM   1262  OE1 GLU A 245      36.941  25.586  57.308  1.00 27.58           O  
ATOM   1263  OE2 GLU A 245      35.300  26.573  58.388  1.00 26.67           O  
ATOM   1264  N   ASP A 246      39.601  21.989  58.085  1.00 25.34           N  
ATOM   1265  CA  ASP A 246      40.268  20.703  57.927  1.00 26.31           C  
ATOM   1266  C   ASP A 246      41.783  20.865  57.972  1.00 26.08           C  
ATOM   1267  O   ASP A 246      42.451  20.164  58.729  1.00 26.16           O  
ATOM   1268  CB  ASP A 246      39.853  20.021  56.613  1.00 26.95           C  
ATOM   1269  CG  ASP A 246      38.362  19.722  56.548  1.00 29.37           C  
ATOM   1270  OD1 ASP A 246      37.689  19.733  57.611  1.00 33.10           O  
ATOM   1271  OD2 ASP A 246      37.861  19.461  55.432  1.00 31.35           O  
ATOM   1272  N   TRP A 247      42.315  21.795  57.172  1.00 25.68           N  
ATOM   1273  CA  TRP A 247      43.747  22.034  57.144  1.00 25.33           C  
ATOM   1274  C   TRP A 247      44.265  22.256  58.557  1.00 25.69           C  
ATOM   1275  O   TRP A 247      45.218  21.601  58.984  1.00 25.35           O  
ATOM   1276  CB  TRP A 247      44.144  23.214  56.244  1.00 24.41           C  
ATOM   1277  CG  TRP A 247      45.633  23.415  56.296  1.00 23.69           C  
ATOM   1278  CD1 TRP A 247      46.321  24.281  57.102  1.00 22.44           C  
ATOM   1279  CD2 TRP A 247      46.625  22.674  55.565  1.00 23.28           C  
ATOM   1280  NE1 TRP A 247      47.676  24.144  56.898  1.00 22.30           N  
ATOM   1281  CE2 TRP A 247      47.893  23.166  55.962  1.00 23.42           C  
ATOM   1282  CE3 TRP A 247      46.567  21.650  54.603  1.00 22.14           C  
ATOM   1283  CZ2 TRP A 247      49.102  22.668  55.420  1.00 23.87           C  
ATOM   1284  CZ3 TRP A 247      47.770  21.155  54.069  1.00 22.57           C  
ATOM   1285  CH2 TRP A 247      49.018  21.666  54.482  1.00 22.46           C  
ATOM   1286  N   ILE A 248      43.618  23.168  59.280  1.00 26.34           N  
ATOM   1287  CA  ILE A 248      44.049  23.503  60.634  1.00 26.78           C  
ATOM   1288  C   ILE A 248      44.038  22.274  61.528  1.00 27.11           C  
ATOM   1289  O   ILE A 248      45.020  22.002  62.206  1.00 27.94           O  
ATOM   1290  CB  ILE A 248      43.216  24.660  61.256  1.00 26.79           C  
ATOM   1291  CG1 ILE A 248      43.533  25.989  60.561  1.00 26.79           C  
ATOM   1292  CG2 ILE A 248      43.457  24.777  62.780  1.00 26.56           C  
ATOM   1293  CD1 ILE A 248      44.964  26.444  60.685  1.00 26.53           C  
ATOM   1294  N   SER A 249      42.942  21.526  61.502  1.00 27.17           N  
ATOM   1295  CA  SER A 249      42.772  20.363  62.355  1.00 27.37           C  
ATOM   1296  C   SER A 249      43.821  19.268  62.094  1.00 27.91           C  
ATOM   1297  O   SER A 249      44.327  18.652  63.027  1.00 27.86           O  
ATOM   1298  CB  SER A 249      41.368  19.803  62.162  1.00 27.34           C  
ATOM   1299  OG  SER A 249      41.102  18.777  63.091  1.00 26.86           O  
ATOM   1300  N   ARG A 250      44.154  19.031  60.830  1.00 28.18           N  
ATOM   1301  CA  ARG A 250      45.063  17.949  60.501  1.00 28.73           C  
ATOM   1302  C   ARG A 250      46.534  18.346  60.391  1.00 28.68           C  
ATOM   1303  O   ARG A 250      47.413  17.543  60.708  1.00 28.54           O  
ATOM   1304  CB  ARG A 250      44.602  17.210  59.242  1.00 29.04           C  
ATOM   1305  CG  ARG A 250      43.499  16.204  59.525  1.00 30.78           C  
ATOM   1306  CD  ARG A 250      43.872  14.849  58.980  1.00 33.45           C  
ATOM   1307  NE  ARG A 250      43.148  14.528  57.754  1.00 35.70           N  
ATOM   1308  CZ  ARG A 250      43.696  13.945  56.692  1.00 36.64           C  
ATOM   1309  NH1 ARG A 250      42.946  13.681  55.627  1.00 37.48           N  
ATOM   1310  NH2 ARG A 250      44.994  13.643  56.685  1.00 36.33           N  
ATOM   1311  N   ARG A 251      46.804  19.572  59.950  1.00 28.69           N  
ATOM   1312  CA  ARG A 251      48.183  19.983  59.672  1.00 28.85           C  
ATOM   1313  C   ARG A 251      48.601  21.351  60.250  1.00 29.18           C  
ATOM   1314  O   ARG A 251      49.757  21.766  60.110  1.00 29.14           O  
ATOM   1315  CB  ARG A 251      48.474  19.879  58.167  1.00 28.61           C  
ATOM   1316  CG  ARG A 251      48.489  18.416  57.643  1.00 28.68           C  
ATOM   1317  CD  ARG A 251      49.393  18.219  56.428  1.00 29.27           C  
ATOM   1318  NE  ARG A 251      50.743  18.726  56.684  1.00 31.06           N  
ATOM   1319  CZ  ARG A 251      51.712  18.829  55.777  1.00 32.23           C  
ATOM   1320  NH1 ARG A 251      51.518  18.447  54.519  1.00 32.93           N  
ATOM   1321  NH2 ARG A 251      52.891  19.319  56.137  1.00 33.53           N  
ATOM   1322  N   GLY A 252      47.669  22.030  60.919  1.00 29.49           N  
ATOM   1323  CA  GLY A 252      47.906  23.372  61.449  1.00 30.14           C  
ATOM   1324  C   GLY A 252      49.009  23.472  62.493  1.00 30.59           C  
ATOM   1325  O   GLY A 252      49.532  24.553  62.743  1.00 30.86           O  
ATOM   1326  N   SER A 253      49.355  22.344  63.106  1.00 30.92           N  
ATOM   1327  CA  SER A 253      50.446  22.279  64.077  1.00 31.12           C  
ATOM   1328  C   SER A 253      51.576  21.389  63.565  1.00 31.43           C  
ATOM   1329  O   SER A 253      52.445  20.974  64.330  1.00 31.21           O  
ATOM   1330  CB  SER A 253      49.937  21.763  65.424  1.00 30.83           C  
ATOM   1331  OG  SER A 253      48.924  22.614  65.930  1.00 30.73           O  
ATOM   1332  N   ASP A 254      51.552  21.113  62.264  1.00 32.00           N  
ATOM   1333  CA  ASP A 254      52.549  20.267  61.625  1.00 32.59           C  
ATOM   1334  C   ASP A 254      53.750  21.117  61.241  1.00 32.33           C  
ATOM   1335  O   ASP A 254      53.679  21.890  60.277  1.00 31.79           O  
ATOM   1336  CB  ASP A 254      51.958  19.603  60.378  1.00 33.10           C  
ATOM   1337  CG  ASP A 254      52.505  18.216  60.134  1.00 34.92           C  
ATOM   1338  OD1 ASP A 254      53.568  17.869  60.706  1.00 36.91           O  
ATOM   1339  OD2 ASP A 254      51.859  17.467  59.364  1.00 36.39           O  
ATOM   1340  N   PRO A 255      54.862  20.978  61.989  1.00 32.38           N  
ATOM   1341  CA  PRO A 255      56.030  21.791  61.667  1.00 32.38           C  
ATOM   1342  C   PRO A 255      56.515  21.430  60.275  1.00 32.39           C  
ATOM   1343  O   PRO A 255      56.592  20.242  59.932  1.00 32.45           O  
ATOM   1344  CB  PRO A 255      57.069  21.369  62.718  1.00 32.16           C  
ATOM   1345  CG  PRO A 255      56.292  20.701  63.797  1.00 32.26           C  
ATOM   1346  CD  PRO A 255      55.129  20.061  63.112  1.00 32.34           C  
ATOM   1347  N   VAL A 256      56.793  22.452  59.476  1.00 32.37           N  
ATOM   1348  CA  VAL A 256      57.312  22.259  58.133  1.00 32.58           C  
ATOM   1349  C   VAL A 256      58.747  22.780  58.065  1.00 32.44           C  
ATOM   1350  O   VAL A 256      59.097  23.766  58.716  1.00 32.11           O  
ATOM   1351  CB  VAL A 256      56.403  22.915  57.066  1.00 32.85           C  
ATOM   1352  CG1 VAL A 256      54.996  22.332  57.147  1.00 33.69           C  
ATOM   1353  CG2 VAL A 256      56.323  24.409  57.268  1.00 33.28           C  
ATOM   1354  N   SER A 257      59.571  22.093  57.285  1.00 32.38           N  
ATOM   1355  CA  SER A 257      60.991  22.378  57.209  1.00 32.24           C  
ATOM   1356  C   SER A 257      61.230  23.732  56.569  1.00 32.44           C  
ATOM   1357  O   SER A 257      60.692  24.031  55.500  1.00 32.33           O  
ATOM   1358  CB  SER A 257      61.691  21.288  56.402  1.00 32.20           C  
ATOM   1359  OG  SER A 257      63.092  21.467  56.416  1.00 32.39           O  
ATOM   1360  N   LEU A 258      62.034  24.557  57.232  1.00 32.74           N  
ATOM   1361  CA  LEU A 258      62.448  25.836  56.664  1.00 33.03           C  
ATOM   1362  C   LEU A 258      63.161  25.582  55.350  1.00 33.00           C  
ATOM   1363  O   LEU A 258      62.933  26.270  54.353  1.00 33.09           O  
ATOM   1364  CB  LEU A 258      63.366  26.593  57.625  1.00 32.96           C  
ATOM   1365  CG  LEU A 258      62.711  27.209  58.862  1.00 33.64           C  
ATOM   1366  CD1 LEU A 258      63.779  27.578  59.890  1.00 34.61           C  
ATOM   1367  CD2 LEU A 258      61.853  28.417  58.503  1.00 33.15           C  
ATOM   1368  N   ASP A 259      64.008  24.562  55.363  1.00 33.16           N  
ATOM   1369  CA  ASP A 259      64.748  24.143  54.195  1.00 33.57           C  
ATOM   1370  C   ASP A 259      63.830  23.852  53.006  1.00 33.54           C  
ATOM   1371  O   ASP A 259      63.977  24.452  51.940  1.00 33.49           O  
ATOM   1372  CB  ASP A 259      65.579  22.924  54.546  1.00 33.69           C  
ATOM   1373  CG  ASP A 259      66.760  22.769  53.652  1.00 34.81           C  
ATOM   1374  OD1 ASP A 259      66.712  21.862  52.792  1.00 36.88           O  
ATOM   1375  OD2 ASP A 259      67.724  23.555  53.803  1.00 34.85           O  
ATOM   1376  N   GLU A 260      62.877  22.947  53.198  1.00 33.76           N  
ATOM   1377  CA  GLU A 260      61.875  22.648  52.174  1.00 34.13           C  
ATOM   1378  C   GLU A 260      61.127  23.898  51.736  1.00 34.16           C  
ATOM   1379  O   GLU A 260      60.843  24.079  50.550  1.00 34.23           O  
ATOM   1380  CB  GLU A 260      60.881  21.599  52.669  1.00 34.20           C  
ATOM   1381  CG  GLU A 260      61.396  20.163  52.586  1.00 35.02           C  
ATOM   1382  CD  GLU A 260      60.276  19.145  52.687  1.00 36.36           C  
ATOM   1383  OE1 GLU A 260      59.536  19.170  53.687  1.00 36.81           O  
ATOM   1384  OE2 GLU A 260      60.128  18.320  51.765  1.00 37.66           O  
ATOM   1385  N   LEU A 261      60.830  24.766  52.698  1.00 34.25           N  
ATOM   1386  CA  LEU A 261      60.105  25.998  52.435  1.00 34.27           C  
ATOM   1387  C   LEU A 261      60.876  26.942  51.494  1.00 34.51           C  
ATOM   1388  O   LEU A 261      60.274  27.557  50.609  1.00 34.34           O  
ATOM   1389  CB  LEU A 261      59.755  26.677  53.763  1.00 34.37           C  
ATOM   1390  CG  LEU A 261      58.367  27.284  54.014  1.00 34.02           C  
ATOM   1391  CD1 LEU A 261      57.232  26.450  53.450  1.00 32.72           C  
ATOM   1392  CD2 LEU A 261      58.180  27.489  55.512  1.00 34.19           C  
ATOM   1393  N   ASN A 262      62.199  27.039  51.673  1.00 34.73           N  
ATOM   1394  CA  ASN A 262      63.049  27.866  50.799  1.00 34.88           C  
ATOM   1395  C   ASN A 262      63.052  27.363  49.358  1.00 34.86           C  
ATOM   1396  O   ASN A 262      63.038  28.154  48.412  1.00 34.90           O  
ATOM   1397  CB  ASN A 262      64.492  27.958  51.323  1.00 34.88           C  
ATOM   1398  CG  ASN A 262      65.355  28.966  50.527  1.00 35.85           C  
ATOM   1399  OD1 ASN A 262      64.911  30.077  50.203  1.00 35.75           O  
ATOM   1400  ND2 ASN A 262      66.592  28.575  50.219  1.00 35.70           N  
ATOM   1401  N   GLN A 263      63.059  26.045  49.197  1.00 34.95           N  
ATOM   1402  CA  GLN A 263      63.053  25.435  47.874  1.00 35.06           C  
ATOM   1403  C   GLN A 263      61.756  25.726  47.123  1.00 35.14           C  
ATOM   1404  O   GLN A 263      61.769  25.870  45.894  1.00 35.36           O  
ATOM   1405  CB  GLN A 263      63.296  23.936  47.971  1.00 35.04           C  
ATOM   1406  CG  GLN A 263      64.720  23.584  48.338  1.00 35.60           C  
ATOM   1407  CD  GLN A 263      64.821  22.222  48.982  1.00 36.58           C  
ATOM   1408  OE1 GLN A 263      64.727  21.191  48.314  1.00 36.22           O  
ATOM   1409  NE2 GLN A 263      65.015  22.209  50.294  1.00 38.05           N  
ATOM   1410  N   LEU A 264      60.647  25.818  47.858  1.00 34.86           N  
ATOM   1411  CA  LEU A 264      59.372  26.178  47.257  1.00 34.49           C  
ATOM   1412  C   LEU A 264      59.416  27.640  46.871  1.00 34.40           C  
ATOM   1413  O   LEU A 264      59.110  27.993  45.740  1.00 34.41           O  
ATOM   1414  CB  LEU A 264      58.213  25.950  48.228  1.00 34.47           C  
ATOM   1415  CG  LEU A 264      56.797  25.626  47.711  1.00 34.69           C  
ATOM   1416  CD1 LEU A 264      55.740  26.330  48.563  1.00 32.91           C  
ATOM   1417  CD2 LEU A 264      56.564  25.918  46.218  1.00 34.01           C  
ATOM   1418  N   LEU A 265      59.810  28.481  47.823  1.00 34.36           N  
ATOM   1419  CA  LEU A 265      59.748  29.928  47.667  1.00 34.20           C  
ATOM   1420  C   LEU A 265      61.021  30.467  47.026  1.00 34.44           C  
ATOM   1421  O   LEU A 265      61.821  31.173  47.656  1.00 34.80           O  
ATOM   1422  CB  LEU A 265      59.480  30.598  49.016  1.00 34.11           C  
ATOM   1423  CG  LEU A 265      58.221  30.192  49.796  1.00 34.14           C  
ATOM   1424  CD1 LEU A 265      57.957  31.187  50.917  1.00 34.10           C  
ATOM   1425  CD2 LEU A 265      56.989  30.039  48.902  1.00 33.42           C  
ATOM   1426  N   SER A 266      61.186  30.132  45.755  1.00 34.30           N  
ATOM   1427  CA  SER A 266      62.405  30.403  45.027  1.00 34.26           C  
ATOM   1428  C   SER A 266      62.065  30.492  43.551  1.00 34.38           C  
ATOM   1429  O   SER A 266      61.263  29.702  43.053  1.00 34.46           O  
ATOM   1430  CB  SER A 266      63.404  29.272  45.278  1.00 34.24           C  
ATOM   1431  OG  SER A 266      64.197  29.007  44.139  1.00 33.83           O  
ATOM   1432  N   LYS A 267      62.671  31.456  42.859  1.00 34.38           N  
ATOM   1433  CA  LYS A 267      62.472  31.634  41.417  1.00 34.30           C  
ATOM   1434  C   LYS A 267      62.828  30.341  40.700  1.00 34.11           C  
ATOM   1435  O   LYS A 267      62.241  30.000  39.674  1.00 34.14           O  
ATOM   1436  CB  LYS A 267      63.357  32.765  40.885  1.00 34.33           C  
ATOM   1437  CG  LYS A 267      63.396  34.021  41.751  1.00 34.71           C  
ATOM   1438  CD  LYS A 267      62.313  35.015  41.364  1.00 35.26           C  
ATOM   1439  CE  LYS A 267      62.223  36.157  42.368  1.00 35.56           C  
ATOM   1440  NZ  LYS A 267      61.536  35.752  43.632  1.00 35.24           N  
ATOM   1441  N   ASP A 268      63.791  29.629  41.279  1.00 33.96           N  
ATOM   1442  CA  ASP A 268      64.342  28.390  40.745  1.00 33.90           C  
ATOM   1443  C   ASP A 268      63.315  27.260  40.685  1.00 33.53           C  
ATOM   1444  O   ASP A 268      63.390  26.399  39.805  1.00 33.48           O  
ATOM   1445  CB  ASP A 268      65.507  27.952  41.633  1.00 34.11           C  
ATOM   1446  CG  ASP A 268      66.699  27.468  40.847  1.00 34.95           C  
ATOM   1447  OD1 ASP A 268      67.674  27.028  41.500  1.00 35.93           O  
ATOM   1448  OD2 ASP A 268      66.675  27.527  39.594  1.00 35.78           O  
ATOM   1449  N   PHE A 269      62.365  27.277  41.625  1.00 33.01           N  
ATOM   1450  CA  PHE A 269      61.340  26.237  41.765  1.00 32.49           C  
ATOM   1451  C   PHE A 269      60.593  25.910  40.463  1.00 32.67           C  
ATOM   1452  O   PHE A 269      60.220  24.758  40.225  1.00 32.53           O  
ATOM   1453  CB  PHE A 269      60.340  26.608  42.875  1.00 32.03           C  
ATOM   1454  CG  PHE A 269      59.110  25.746  42.896  1.00 30.71           C  
ATOM   1455  CD1 PHE A 269      59.115  24.526  43.558  1.00 29.30           C  
ATOM   1456  CD2 PHE A 269      57.957  26.140  42.224  1.00 29.98           C  
ATOM   1457  CE1 PHE A 269      57.987  23.714  43.567  1.00 28.97           C  
ATOM   1458  CE2 PHE A 269      56.824  25.330  42.220  1.00 30.31           C  
ATOM   1459  CZ  PHE A 269      56.841  24.112  42.900  1.00 29.86           C  
ATOM   1460  N   SER A 270      60.377  26.925  39.633  1.00 32.90           N  
ATOM   1461  CA  SER A 270      59.672  26.761  38.361  1.00 33.39           C  
ATOM   1462  C   SER A 270      60.301  25.732  37.397  1.00 33.67           C  
ATOM   1463  O   SER A 270      59.584  25.051  36.655  1.00 33.64           O  
ATOM   1464  CB  SER A 270      59.534  28.110  37.661  1.00 33.28           C  
ATOM   1465  OG  SER A 270      58.688  27.990  36.540  1.00 33.77           O  
ATOM   1466  N   LYS A 271      61.631  25.627  37.414  1.00 33.90           N  
ATOM   1467  CA  LYS A 271      62.353  24.736  36.505  1.00 34.17           C  
ATOM   1468  C   LYS A 271      62.825  23.437  37.171  1.00 34.23           C  
ATOM   1469  O   LYS A 271      63.696  22.749  36.644  1.00 34.13           O  
ATOM   1470  CB  LYS A 271      63.533  25.471  35.850  1.00 34.28           C  
ATOM   1471  CG  LYS A 271      63.127  26.526  34.814  1.00 34.53           C  
ATOM   1472  CD  LYS A 271      64.310  26.998  33.951  1.00 34.53           C  
ATOM   1473  CE  LYS A 271      65.153  28.074  34.632  1.00 35.25           C  
ATOM   1474  NZ  LYS A 271      64.345  29.253  35.071  1.00 35.35           N  
HETATM 1475  N   MSE A 272      62.245  23.102  38.322  1.00 34.46           N  
HETATM 1476  CA  MSE A 272      62.531  21.834  38.993  1.00 34.65           C  
HETATM 1477  C   MSE A 272      61.788  20.672  38.323  1.00 33.92           C  
HETATM 1478  O   MSE A 272      60.699  20.864  37.783  1.00 33.57           O  
HETATM 1479  CB  MSE A 272      62.168  21.910  40.481  1.00 35.72           C  
HETATM 1480  CG  MSE A 272      63.088  22.792  41.332  1.00 37.72           C  
HETATM 1481 SE   MSE A 272      64.958  22.217  41.335  1.00 45.38          SE  
HETATM 1482  CE  MSE A 272      65.592  23.256  42.869  1.00 40.61           C  
ATOM   1483  N   PRO A 273      62.387  19.463  38.329  1.00 33.42           N  
ATOM   1484  CA  PRO A 273      61.726  18.288  37.765  1.00 33.04           C  
ATOM   1485  C   PRO A 273      60.352  18.017  38.371  1.00 32.79           C  
ATOM   1486  O   PRO A 273      60.131  18.287  39.546  1.00 32.62           O  
ATOM   1487  CB  PRO A 273      62.696  17.156  38.094  1.00 32.89           C  
ATOM   1488  CG  PRO A 273      64.020  17.823  38.154  1.00 32.88           C  
ATOM   1489  CD  PRO A 273      63.738  19.131  38.816  1.00 33.37           C  
ATOM   1490  N   ASP A 274      59.452  17.476  37.552  1.00 32.70           N  
ATOM   1491  CA  ASP A 274      58.047  17.270  37.909  1.00 32.54           C  
ATOM   1492  C   ASP A 274      57.856  16.624  39.276  1.00 32.26           C  
ATOM   1493  O   ASP A 274      56.984  17.042  40.046  1.00 32.52           O  
ATOM   1494  CB  ASP A 274      57.332  16.431  36.837  1.00 32.76           C  
ATOM   1495  CG  ASP A 274      57.255  17.136  35.480  1.00 33.45           C  
ATOM   1496  OD1 ASP A 274      57.530  18.360  35.401  1.00 33.99           O  
ATOM   1497  OD2 ASP A 274      56.912  16.455  34.487  1.00 33.72           O  
ATOM   1498  N   ASP A 275      58.675  15.619  39.578  1.00 31.63           N  
ATOM   1499  CA  ASP A 275      58.558  14.899  40.838  1.00 31.07           C  
ATOM   1500  C   ASP A 275      59.079  15.684  42.039  1.00 30.58           C  
ATOM   1501  O   ASP A 275      58.555  15.547  43.141  1.00 30.75           O  
ATOM   1502  CB  ASP A 275      59.189  13.505  40.746  1.00 31.17           C  
ATOM   1503  CG  ASP A 275      58.248  12.476  40.109  1.00 31.90           C  
ATOM   1504  OD1 ASP A 275      58.539  11.263  40.176  1.00 32.40           O  
ATOM   1505  OD2 ASP A 275      57.207  12.874  39.540  1.00 32.83           O  
ATOM   1506  N   VAL A 276      60.086  16.523  41.820  1.00 29.85           N  
ATOM   1507  CA  VAL A 276      60.579  17.412  42.865  1.00 29.06           C  
ATOM   1508  C   VAL A 276      59.488  18.411  43.226  1.00 28.89           C  
ATOM   1509  O   VAL A 276      59.166  18.586  44.395  1.00 28.68           O  
ATOM   1510  CB  VAL A 276      61.873  18.142  42.440  1.00 28.92           C  
ATOM   1511  CG1 VAL A 276      62.322  19.119  43.513  1.00 28.24           C  
ATOM   1512  CG2 VAL A 276      62.971  17.139  42.152  1.00 28.37           C  
ATOM   1513  N   GLN A 277      58.908  19.039  42.210  1.00 29.08           N  
ATOM   1514  CA  GLN A 277      57.827  20.012  42.393  1.00 29.47           C  
ATOM   1515  C   GLN A 277      56.639  19.375  43.116  1.00 29.71           C  
ATOM   1516  O   GLN A 277      56.241  19.809  44.199  1.00 29.80           O  
ATOM   1517  CB  GLN A 277      57.389  20.587  41.042  1.00 29.33           C  
ATOM   1518  CG  GLN A 277      58.334  21.635  40.455  1.00 29.49           C  
ATOM   1519  CD  GLN A 277      57.943  22.077  39.047  1.00 29.46           C  
ATOM   1520  OE1 GLN A 277      57.173  21.405  38.360  1.00 30.14           O  
ATOM   1521  NE2 GLN A 277      58.482  23.210  38.612  1.00 29.15           N  
ATOM   1522  N   THR A 278      56.101  18.323  42.513  1.00 29.88           N  
ATOM   1523  CA  THR A 278      55.032  17.541  43.104  1.00 30.08           C  
ATOM   1524  C   THR A 278      55.306  17.204  44.572  1.00 30.09           C  
ATOM   1525  O   THR A 278      54.437  17.402  45.427  1.00 30.25           O  
ATOM   1526  CB  THR A 278      54.824  16.258  42.288  1.00 30.12           C  
ATOM   1527  OG1 THR A 278      54.460  16.622  40.948  1.00 30.94           O  
ATOM   1528  CG2 THR A 278      53.734  15.395  42.888  1.00 30.13           C  
ATOM   1529  N   ARG A 279      56.514  16.715  44.860  1.00 29.82           N  
ATOM   1530  CA  ARG A 279      56.862  16.263  46.205  1.00 29.41           C  
ATOM   1531  C   ARG A 279      56.821  17.399  47.225  1.00 29.01           C  
ATOM   1532  O   ARG A 279      56.282  17.233  48.315  1.00 29.10           O  
ATOM   1533  CB  ARG A 279      58.228  15.581  46.211  1.00 29.40           C  
ATOM   1534  CG  ARG A 279      58.557  14.834  47.499  1.00 30.73           C  
ATOM   1535  CD  ARG A 279      57.536  13.735  47.834  1.00 32.52           C  
ATOM   1536  NE  ARG A 279      57.647  12.571  46.952  1.00 33.30           N  
ATOM   1537  CZ  ARG A 279      56.965  11.438  47.113  1.00 33.92           C  
ATOM   1538  NH1 ARG A 279      57.134  10.435  46.264  1.00 33.56           N  
ATOM   1539  NH2 ARG A 279      56.111  11.302  48.123  1.00 34.35           N  
ATOM   1540  N   LEU A 280      57.374  18.548  46.844  1.00 28.22           N  
ATOM   1541  CA  LEU A 280      57.400  19.732  47.678  1.00 27.45           C  
ATOM   1542  C   LEU A 280      56.003  20.269  47.958  1.00 27.36           C  
ATOM   1543  O   LEU A 280      55.672  20.606  49.098  1.00 27.20           O  
ATOM   1544  CB  LEU A 280      58.253  20.820  47.023  1.00 27.17           C  
ATOM   1545  CG  LEU A 280      59.766  20.783  47.242  1.00 26.67           C  
ATOM   1546  CD1 LEU A 280      60.447  21.708  46.242  1.00 26.02           C  
ATOM   1547  CD2 LEU A 280      60.126  21.170  48.679  1.00 24.96           C  
ATOM   1548  N   LEU A 281      55.187  20.350  46.915  1.00 27.15           N  
ATOM   1549  CA  LEU A 281      53.826  20.842  47.070  1.00 27.07           C  
ATOM   1550  C   LEU A 281      53.031  19.965  48.039  1.00 27.06           C  
ATOM   1551  O   LEU A 281      52.415  20.466  48.975  1.00 26.52           O  
ATOM   1552  CB  LEU A 281      53.129  20.940  45.714  1.00 26.78           C  
ATOM   1553  CG  LEU A 281      53.652  22.032  44.777  1.00 26.82           C  
ATOM   1554  CD1 LEU A 281      52.935  21.968  43.445  1.00 27.80           C  
ATOM   1555  CD2 LEU A 281      53.497  23.427  45.379  1.00 27.60           C  
ATOM   1556  N   ALA A 282      53.075  18.654  47.812  1.00 27.11           N  
ATOM   1557  CA  ALA A 282      52.356  17.699  48.643  1.00 27.18           C  
ATOM   1558  C   ALA A 282      52.796  17.826  50.087  1.00 27.03           C  
ATOM   1559  O   ALA A 282      51.970  17.918  50.987  1.00 26.96           O  
ATOM   1560  CB  ALA A 282      52.601  16.280  48.149  1.00 27.23           C  
ATOM   1561  N   SER A 283      54.113  17.840  50.274  1.00 27.06           N  
ATOM   1562  CA  SER A 283      54.763  17.852  51.575  1.00 27.07           C  
ATOM   1563  C   SER A 283      54.435  19.106  52.385  1.00 27.20           C  
ATOM   1564  O   SER A 283      54.127  19.025  53.568  1.00 27.33           O  
ATOM   1565  CB  SER A 283      56.271  17.719  51.374  1.00 26.97           C  
ATOM   1566  OG  SER A 283      56.985  18.110  52.523  1.00 28.08           O  
ATOM   1567  N   ILE A 284      54.487  20.260  51.732  1.00 27.33           N  
ATOM   1568  CA  ILE A 284      54.244  21.532  52.380  1.00 27.37           C  
ATOM   1569  C   ILE A 284      52.747  21.879  52.437  1.00 27.66           C  
ATOM   1570  O   ILE A 284      52.258  22.373  53.466  1.00 28.17           O  
ATOM   1571  CB  ILE A 284      55.073  22.672  51.688  1.00 27.60           C  
ATOM   1572  CG1 ILE A 284      56.586  22.469  51.918  1.00 27.65           C  
ATOM   1573  CG2 ILE A 284      54.664  24.051  52.185  1.00 27.01           C  
ATOM   1574  CD1 ILE A 284      57.485  23.284  50.992  1.00 27.11           C  
ATOM   1575  N   LEU A 285      52.016  21.626  51.354  1.00 27.26           N  
ATOM   1576  CA  LEU A 285      50.669  22.189  51.233  1.00 27.07           C  
ATOM   1577  C   LEU A 285      49.503  21.220  51.042  1.00 27.18           C  
ATOM   1578  O   LEU A 285      48.369  21.663  50.890  1.00 27.18           O  
ATOM   1579  CB  LEU A 285      50.625  23.244  50.125  1.00 26.76           C  
ATOM   1580  CG  LEU A 285      51.552  24.457  50.216  1.00 26.71           C  
ATOM   1581  CD1 LEU A 285      51.452  25.257  48.929  1.00 25.12           C  
ATOM   1582  CD2 LEU A 285      51.274  25.344  51.434  1.00 25.56           C  
ATOM   1583  N   GLN A 286      49.756  19.919  51.045  1.00 27.44           N  
ATOM   1584  CA  GLN A 286      48.666  18.969  50.847  1.00 27.95           C  
ATOM   1585  C   GLN A 286      48.362  18.240  52.154  1.00 28.50           C  
ATOM   1586  O   GLN A 286      49.279  17.798  52.850  1.00 28.34           O  
ATOM   1587  CB  GLN A 286      48.983  18.000  49.708  1.00 27.66           C  
ATOM   1588  CG  GLN A 286      47.755  17.494  48.965  1.00 27.37           C  
ATOM   1589  CD  GLN A 286      47.113  16.289  49.624  1.00 26.71           C  
ATOM   1590  OE1 GLN A 286      45.894  16.136  49.599  1.00 27.05           O  
ATOM   1591  NE2 GLN A 286      47.931  15.428  50.212  1.00 25.25           N  
ATOM   1592  N   ILE A 287      47.072  18.129  52.476  1.00 29.27           N  
ATOM   1593  CA  ILE A 287      46.621  17.656  53.790  1.00 30.16           C  
ATOM   1594  C   ILE A 287      47.058  16.228  54.121  1.00 31.20           C  
ATOM   1595  O   ILE A 287      47.143  15.846  55.286  1.00 31.45           O  
ATOM   1596  CB  ILE A 287      45.077  17.813  53.965  1.00 30.34           C  
ATOM   1597  CG1 ILE A 287      44.725  17.947  55.456  1.00 30.50           C  
ATOM   1598  CG2 ILE A 287      44.285  16.685  53.231  1.00 28.72           C  
ATOM   1599  CD1 ILE A 287      43.291  18.374  55.748  1.00 30.48           C  
ATOM   1600  N   ASP A 288      47.342  15.443  53.090  1.00 32.12           N  
ATOM   1601  CA  ASP A 288      47.701  14.054  53.279  1.00 32.85           C  
ATOM   1602  C   ASP A 288      49.142  13.765  52.854  1.00 33.10           C  
ATOM   1603  O   ASP A 288      49.539  12.603  52.820  1.00 33.28           O  
ATOM   1604  CB  ASP A 288      46.702  13.150  52.538  1.00 33.24           C  
ATOM   1605  CG  ASP A 288      45.472  12.814  53.382  1.00 34.86           C  
ATOM   1606  OD1 ASP A 288      45.646  12.169  54.442  1.00 37.61           O  
ATOM   1607  OD2 ASP A 288      44.336  13.178  52.991  1.00 34.69           O  
ATOM   1608  N   LYS A 289      49.888  14.798  52.519  1.00 20.00           N  
ATOM   1609  CA  LYS A 289      51.238  14.667  52.031  1.00 20.00           C  
ATOM   1610  C   LYS A 289      51.402  13.807  50.814  1.00 20.00           C  
ATOM   1611  O   LYS A 289      52.352  13.104  50.654  1.00 33.73           O  
ATOM   1612  CB  LYS A 289      52.204  14.283  53.117  1.00 20.00           C  
ATOM   1613  CG  LYS A 289      52.170  15.212  54.261  1.00 20.00           C  
ATOM   1614  CD  LYS A 289      53.502  15.551  54.716  1.00 20.00           C  
ATOM   1615  CE  LYS A 289      53.540  15.504  56.211  1.00 20.00           C  
ATOM   1616  NZ  LYS A 289      54.813  16.008  56.783  1.00 20.00           N  
ATOM   1617  N   ASP A 290      50.461  13.937  49.924  1.00 33.15           N  
ATOM   1618  CA  ASP A 290      50.231  12.976  48.852  1.00 32.78           C  
ATOM   1619  C   ASP A 290      50.473  13.553  47.450  1.00 32.74           C  
ATOM   1620  O   ASP A 290      49.661  14.333  46.945  1.00 32.52           O  
ATOM   1621  CB  ASP A 290      48.808  12.404  48.977  1.00 32.69           C  
ATOM   1622  CG  ASP A 290      48.508  11.327  47.954  1.00 31.71           C  
ATOM   1623  OD1 ASP A 290      49.448  10.827  47.304  1.00 29.80           O  
ATOM   1624  OD2 ASP A 290      47.317  10.981  47.812  1.00 31.17           O  
ATOM   1625  N   PRO A 291      51.597  13.158  46.820  1.00 32.76           N  
ATOM   1626  CA  PRO A 291      51.983  13.543  45.460  1.00 32.62           C  
ATOM   1627  C   PRO A 291      50.874  13.378  44.426  1.00 32.54           C  
ATOM   1628  O   PRO A 291      50.801  14.152  43.478  1.00 32.50           O  
ATOM   1629  CB  PRO A 291      53.113  12.567  45.135  1.00 32.54           C  
ATOM   1630  CG  PRO A 291      53.711  12.261  46.428  1.00 32.85           C  
ATOM   1631  CD  PRO A 291      52.609  12.285  47.441  1.00 32.76           C  
ATOM   1632  N   HIS A 292      50.023  12.375  44.613  1.00 32.58           N  
ATOM   1633  CA  HIS A 292      49.006  12.024  43.623  1.00 32.71           C  
ATOM   1634  C   HIS A 292      47.831  13.000  43.586  1.00 32.48           C  
ATOM   1635  O   HIS A 292      47.076  13.028  42.612  1.00 32.48           O  
ATOM   1636  CB  HIS A 292      48.514  10.586  43.833  1.00 32.79           C  
ATOM   1637  CG  HIS A 292      49.492   9.545  43.383  1.00 33.86           C  
ATOM   1638  ND1 HIS A 292      50.549   9.127  44.166  1.00 34.34           N  
ATOM   1639  CD2 HIS A 292      49.572   8.835  42.231  1.00 34.49           C  
ATOM   1640  CE1 HIS A 292      51.238   8.206  43.515  1.00 34.58           C  
ATOM   1641  NE2 HIS A 292      50.666   8.009  42.339  1.00 34.44           N  
ATOM   1642  N   LYS A 293      47.683  13.798  44.640  1.00 32.11           N  
ATOM   1643  CA  LYS A 293      46.636  14.818  44.688  1.00 31.84           C  
ATOM   1644  C   LYS A 293      47.118  16.161  44.121  1.00 31.21           C  
ATOM   1645  O   LYS A 293      46.351  17.123  44.050  1.00 31.35           O  
ATOM   1646  CB  LYS A 293      46.104  14.985  46.120  1.00 32.06           C  
ATOM   1647  CG  LYS A 293      45.625  13.685  46.787  1.00 33.67           C  
ATOM   1648  CD  LYS A 293      44.254  13.221  46.273  1.00 36.48           C  
ATOM   1649  CE  LYS A 293      43.128  14.049  46.903  1.00 39.16           C  
ATOM   1650  NZ  LYS A 293      41.815  13.920  46.194  1.00 40.38           N  
ATOM   1651  N   VAL A 294      48.381  16.213  43.705  1.00 30.41           N  
ATOM   1652  CA  VAL A 294      48.973  17.426  43.137  1.00 29.79           C  
ATOM   1653  C   VAL A 294      49.081  17.316  41.616  1.00 29.50           C  
ATOM   1654  O   VAL A 294      49.315  16.238  41.082  1.00 29.15           O  
ATOM   1655  CB  VAL A 294      50.353  17.744  43.782  1.00 29.81           C  
ATOM   1656  CG1 VAL A 294      51.133  18.763  42.963  1.00 29.78           C  
ATOM   1657  CG2 VAL A 294      50.171  18.251  45.205  1.00 29.11           C  
ATOM   1658  N   ASP A 295      48.896  18.437  40.923  1.00 29.37           N  
ATOM   1659  CA  ASP A 295      48.925  18.450  39.461  1.00 29.23           C  
ATOM   1660  C   ASP A 295      49.764  19.619  38.944  1.00 29.01           C  
ATOM   1661  O   ASP A 295      49.240  20.690  38.656  1.00 28.87           O  
ATOM   1662  CB  ASP A 295      47.495  18.498  38.910  1.00 29.23           C  
ATOM   1663  CG  ASP A 295      47.424  18.213  37.419  1.00 30.05           C  
ATOM   1664  OD1 ASP A 295      48.424  17.724  36.842  1.00 31.15           O  
ATOM   1665  OD2 ASP A 295      46.355  18.474  36.818  1.00 30.44           O  
ATOM   1666  N   ILE A 296      51.073  19.400  38.829  1.00 28.93           N  
ATOM   1667  CA  ILE A 296      52.012  20.457  38.435  1.00 28.95           C  
ATOM   1668  C   ILE A 296      51.803  20.958  37.008  1.00 29.00           C  
ATOM   1669  O   ILE A 296      52.398  21.955  36.597  1.00 28.83           O  
ATOM   1670  CB  ILE A 296      53.491  20.039  38.628  1.00 28.95           C  
ATOM   1671  CG1 ILE A 296      53.732  18.602  38.152  1.00 28.93           C  
ATOM   1672  CG2 ILE A 296      53.905  20.211  40.073  1.00 28.96           C  
ATOM   1673  CD1 ILE A 296      54.161  18.505  36.714  1.00 28.90           C  
ATOM   1674  N   LYS A 297      50.957  20.258  36.260  1.00 29.00           N  
ATOM   1675  CA  LYS A 297      50.560  20.710  34.937  1.00 29.32           C  
ATOM   1676  C   LYS A 297      49.651  21.943  35.022  1.00 29.28           C  
ATOM   1677  O   LYS A 297      49.432  22.627  34.023  1.00 29.28           O  
ATOM   1678  CB  LYS A 297      49.911  19.564  34.137  1.00 29.42           C  
ATOM   1679  CG  LYS A 297      50.798  18.936  33.044  1.00 29.84           C  
ATOM   1680  CD  LYS A 297      52.254  18.696  33.482  1.00 31.26           C  
ATOM   1681  CE  LYS A 297      53.210  18.635  32.277  1.00 32.55           C  
ATOM   1682  NZ  LYS A 297      54.647  18.791  32.667  1.00 32.64           N  
ATOM   1683  N   LYS A 298      49.144  22.226  36.221  1.00 29.43           N  
ATOM   1684  CA  LYS A 298      48.320  23.415  36.463  1.00 29.57           C  
ATOM   1685  C   LYS A 298      49.123  24.617  36.979  1.00 29.65           C  
ATOM   1686  O   LYS A 298      48.561  25.698  37.158  1.00 29.56           O  
ATOM   1687  CB  LYS A 298      47.161  23.102  37.421  1.00 29.47           C  
ATOM   1688  CG  LYS A 298      46.254  21.954  36.990  1.00 29.75           C  
ATOM   1689  CD  LYS A 298      45.315  22.327  35.849  1.00 30.74           C  
ATOM   1690  CE  LYS A 298      43.967  22.814  36.366  1.00 31.28           C  
ATOM   1691  NZ  LYS A 298      42.919  22.789  35.299  1.00 31.91           N  
ATOM   1692  N   LEU A 299      50.424  24.425  37.217  1.00 29.89           N  
ATOM   1693  CA  LEU A 299      51.330  25.515  37.607  1.00 30.26           C  
ATOM   1694  C   LEU A 299      51.440  26.573  36.516  1.00 30.59           C  
ATOM   1695  O   LEU A 299      51.499  26.240  35.329  1.00 30.76           O  
ATOM   1696  CB  LEU A 299      52.739  24.982  37.908  1.00 30.07           C  
ATOM   1697  CG  LEU A 299      53.054  24.160  39.162  1.00 30.24           C  
ATOM   1698  CD1 LEU A 299      54.548  23.845  39.207  1.00 30.02           C  
ATOM   1699  CD2 LEU A 299      52.631  24.873  40.444  1.00 29.25           C  
ATOM   1700  N   HIS A 300      51.461  27.844  36.913  1.00 31.21           N  
ATOM   1701  CA  HIS A 300      51.783  28.924  35.976  1.00 31.76           C  
ATOM   1702  C   HIS A 300      53.281  29.191  36.089  1.00 32.00           C  
ATOM   1703  O   HIS A 300      53.728  29.950  36.947  1.00 32.17           O  
ATOM   1704  CB  HIS A 300      50.953  30.193  36.224  1.00 31.51           C  
ATOM   1705  CG  HIS A 300      51.439  31.387  35.456  1.00 32.78           C  
ATOM   1706  ND1 HIS A 300      51.328  31.490  34.084  1.00 33.19           N  
ATOM   1707  CD2 HIS A 300      52.062  32.520  35.866  1.00 33.13           C  
ATOM   1708  CE1 HIS A 300      51.856  32.633  33.683  1.00 32.50           C  
ATOM   1709  NE2 HIS A 300      52.309  33.277  34.744  1.00 32.72           N  
ATOM   1710  N   LEU A 301      54.049  28.546  35.216  1.00 32.49           N  
ATOM   1711  CA  LEU A 301      55.508  28.475  35.356  1.00 32.85           C  
ATOM   1712  C   LEU A 301      56.211  29.837  35.268  1.00 33.17           C  
ATOM   1713  O   LEU A 301      57.333  29.995  35.744  1.00 33.13           O  
ATOM   1714  CB  LEU A 301      56.094  27.478  34.343  1.00 32.83           C  
ATOM   1715  CG  LEU A 301      55.559  26.036  34.366  1.00 32.69           C  
ATOM   1716  CD1 LEU A 301      55.854  25.328  33.053  1.00 32.82           C  
ATOM   1717  CD2 LEU A 301      56.105  25.237  35.541  1.00 32.12           C  
ATOM   1718  N   ASP A 302      55.543  30.820  34.673  1.00 33.62           N  
ATOM   1719  CA  ASP A 302      56.101  32.162  34.578  1.00 34.00           C  
ATOM   1720  C   ASP A 302      55.597  33.028  35.734  1.00 34.19           C  
ATOM   1721  O   ASP A 302      55.543  34.257  35.639  1.00 34.21           O  
ATOM   1722  CB  ASP A 302      55.768  32.789  33.222  1.00 33.92           C  
ATOM   1723  CG  ASP A 302      56.995  33.345  32.527  1.00 34.44           C  
ATOM   1724  OD1 ASP A 302      56.878  33.752  31.348  1.00 34.87           O  
ATOM   1725  OD2 ASP A 302      58.082  33.365  33.154  1.00 34.42           O  
ATOM   1726  N   GLY A 303      55.247  32.364  36.832  1.00 34.40           N  
ATOM   1727  CA  GLY A 303      54.662  33.021  37.992  1.00 34.66           C  
ATOM   1728  C   GLY A 303      55.617  33.897  38.773  1.00 34.72           C  
ATOM   1729  O   GLY A 303      56.831  33.683  38.766  1.00 34.82           O  
ATOM   1730  N   LYS A 304      55.043  34.887  39.444  1.00 34.60           N  
ATOM   1731  CA  LYS A 304      55.786  35.806  40.274  1.00 34.81           C  
ATOM   1732  C   LYS A 304      55.583  35.419  41.749  1.00 34.42           C  
ATOM   1733  O   LYS A 304      54.448  35.258  42.206  1.00 34.53           O  
ATOM   1734  CB  LYS A 304      55.280  37.226  40.002  1.00 34.89           C  
ATOM   1735  CG  LYS A 304      56.345  38.314  40.046  1.00 35.71           C  
ATOM   1736  CD  LYS A 304      55.730  39.705  39.846  1.00 35.70           C  
ATOM   1737  CE  LYS A 304      55.788  40.159  38.389  1.00 37.04           C  
ATOM   1738  NZ  LYS A 304      56.985  40.992  38.097  1.00 37.81           N  
ATOM   1739  N   LEU A 305      56.673  35.251  42.491  1.00 33.98           N  
ATOM   1740  CA  LEU A 305      56.569  35.011  43.930  1.00 33.55           C  
ATOM   1741  C   LEU A 305      56.814  36.278  44.747  1.00 33.25           C  
ATOM   1742  O   LEU A 305      57.869  36.895  44.632  1.00 33.22           O  
ATOM   1743  CB  LEU A 305      57.542  33.921  44.373  1.00 33.85           C  
ATOM   1744  CG  LEU A 305      57.663  33.739  45.896  1.00 34.23           C  
ATOM   1745  CD1 LEU A 305      56.576  32.823  46.439  1.00 33.46           C  
ATOM   1746  CD2 LEU A 305      59.035  33.228  46.263  1.00 35.11           C  
ATOM   1747  N   ARG A 306      55.835  36.652  45.570  1.00 32.98           N  
ATOM   1748  CA  ARG A 306      55.940  37.808  46.464  1.00 33.02           C  
ATOM   1749  C   ARG A 306      55.466  37.404  47.857  1.00 32.26           C  
ATOM   1750  O   ARG A 306      54.549  36.602  47.978  1.00 32.45           O  
ATOM   1751  CB  ARG A 306      55.096  38.976  45.946  1.00 32.93           C  
ATOM   1752  CG  ARG A 306      54.834  38.948  44.429  1.00 34.79           C  
ATOM   1753  CD  ARG A 306      54.351  40.282  43.848  1.00 34.49           C  
ATOM   1754  NE  ARG A 306      53.235  40.860  44.598  1.00 37.99           N  
ATOM   1755  CZ  ARG A 306      53.267  42.053  45.195  1.00 40.01           C  
ATOM   1756  NH1 ARG A 306      54.353  42.820  45.118  1.00 40.45           N  
ATOM   1757  NH2 ARG A 306      52.201  42.493  45.859  1.00 41.17           N  
ATOM   1758  N   PHE A 307      56.088  37.952  48.902  1.00 31.58           N  
ATOM   1759  CA  PHE A 307      55.677  37.688  50.293  1.00 30.75           C  
ATOM   1760  C   PHE A 307      55.944  38.873  51.219  1.00 30.55           C  
ATOM   1761  O   PHE A 307      56.803  39.711  50.950  1.00 30.30           O  
ATOM   1762  CB  PHE A 307      56.335  36.412  50.856  1.00 30.81           C  
ATOM   1763  CG  PHE A 307      57.849  36.437  50.842  1.00 30.46           C  
ATOM   1764  CD1 PHE A 307      58.561  36.973  51.910  1.00 30.80           C  
ATOM   1765  CD2 PHE A 307      58.556  35.929  49.759  1.00 30.42           C  
ATOM   1766  CE1 PHE A 307      59.959  37.009  51.900  1.00 31.77           C  
ATOM   1767  CE2 PHE A 307      59.953  35.958  49.737  1.00 31.68           C  
ATOM   1768  CZ  PHE A 307      60.655  36.499  50.813  1.00 31.54           C  
ATOM   1769  N   ALA A 308      55.199  38.940  52.315  1.00 30.30           N  
ATOM   1770  CA  ALA A 308      55.403  39.990  53.308  1.00 30.06           C  
ATOM   1771  C   ALA A 308      55.286  39.386  54.687  1.00 29.81           C  
ATOM   1772  O   ALA A 308      54.357  38.628  54.950  1.00 30.21           O  
ATOM   1773  CB  ALA A 308      54.380  41.087  53.141  1.00 30.05           C  
ATOM   1774  N   SER A 309      56.231  39.722  55.558  1.00 29.37           N  
ATOM   1775  CA  SER A 309      56.186  39.301  56.950  1.00 29.16           C  
ATOM   1776  C   SER A 309      55.759  40.460  57.841  1.00 28.55           C  
ATOM   1777  O   SER A 309      55.865  41.624  57.450  1.00 28.31           O  
ATOM   1778  CB  SER A 309      57.548  38.768  57.393  1.00 29.35           C  
ATOM   1779  OG  SER A 309      58.511  39.806  57.358  1.00 30.01           O  
ATOM   1780  N   HIS A 310      55.275  40.125  59.034  1.00 28.15           N  
ATOM   1781  CA  HIS A 310      54.813  41.102  60.017  1.00 28.15           C  
ATOM   1782  C   HIS A 310      54.941  40.522  61.416  1.00 28.44           C  
ATOM   1783  O   HIS A 310      54.498  39.395  61.661  1.00 28.64           O  
ATOM   1784  CB  HIS A 310      53.348  41.442  59.766  1.00 27.91           C  
ATOM   1785  CG  HIS A 310      52.908  42.716  60.414  1.00 28.23           C  
ATOM   1786  ND1 HIS A 310      52.735  42.838  61.778  1.00 28.22           N  
ATOM   1787  CD2 HIS A 310      52.606  43.926  59.887  1.00 27.21           C  
ATOM   1788  CE1 HIS A 310      52.347  44.068  62.061  1.00 27.12           C  
ATOM   1789  NE2 HIS A 310      52.258  44.747  60.931  1.00 27.40           N  
ATOM   1790  N   GLU A 311      55.544  41.276  62.337  1.00 28.63           N  
ATOM   1791  CA  GLU A 311      55.600  40.839  63.738  1.00 29.14           C  
ATOM   1792  C   GLU A 311      54.186  40.733  64.343  1.00 28.77           C  
ATOM   1793  O   GLU A 311      53.246  41.400  63.875  1.00 28.53           O  
ATOM   1794  CB  GLU A 311      56.540  41.716  64.586  1.00 28.69           C  
ATOM   1795  CG  GLU A 311      55.997  43.086  65.044  1.00 30.02           C  
ATOM   1796  CD  GLU A 311      57.095  43.985  65.654  1.00 31.40           C  
ATOM   1797  OE1 GLU A 311      58.227  43.496  65.897  1.00 35.16           O  
ATOM   1798  OE2 GLU A 311      56.840  45.188  65.896  1.00 34.39           O  
ATOM   1799  N   TYR A 312      54.024  39.866  65.341  1.00 28.40           N  
ATOM   1800  CA  TYR A 312      52.769  39.837  66.099  1.00 28.71           C  
ATOM   1801  C   TYR A 312      52.964  39.488  67.573  1.00 29.01           C  
ATOM   1802  O   TYR A 312      53.904  38.765  67.937  1.00 28.80           O  
ATOM   1803  CB  TYR A 312      51.709  38.935  65.440  1.00 28.28           C  
ATOM   1804  CG  TYR A 312      51.992  37.447  65.508  1.00 27.98           C  
ATOM   1805  CD1 TYR A 312      51.487  36.664  66.548  1.00 26.11           C  
ATOM   1806  CD2 TYR A 312      52.767  36.822  64.530  1.00 28.48           C  
ATOM   1807  CE1 TYR A 312      51.741  35.300  66.615  1.00 26.40           C  
ATOM   1808  CE2 TYR A 312      53.039  35.445  64.596  1.00 29.01           C  
ATOM   1809  CZ  TYR A 312      52.519  34.688  65.640  1.00 27.89           C  
ATOM   1810  OH  TYR A 312      52.778  33.326  65.683  1.00 26.64           O  
ATOM   1811  N   ASP A 313      52.080  40.034  68.410  1.00 29.42           N  
ATOM   1812  CA  ASP A 313      52.037  39.707  69.827  1.00 29.96           C  
ATOM   1813  C   ASP A 313      51.203  38.449  69.995  1.00 30.02           C  
ATOM   1814  O   ASP A 313      50.000  38.432  69.688  1.00 30.07           O  
ATOM   1815  CB  ASP A 313      51.427  40.850  70.638  1.00 30.30           C  
ATOM   1816  CG  ASP A 313      51.435  40.575  72.139  1.00 31.26           C  
ATOM   1817  OD1 ASP A 313      52.128  41.313  72.859  1.00 32.91           O  
ATOM   1818  OD2 ASP A 313      50.763  39.627  72.601  1.00 31.59           O  
ATOM   1819  N   PHE A 314      51.838  37.393  70.478  1.00 29.54           N  
ATOM   1820  CA  PHE A 314      51.122  36.152  70.612  1.00 29.73           C  
ATOM   1821  C   PHE A 314      49.881  36.275  71.520  1.00 29.82           C  
ATOM   1822  O   PHE A 314      48.779  35.841  71.139  1.00 30.01           O  
ATOM   1823  CB  PHE A 314      52.039  35.025  71.074  1.00 29.72           C  
ATOM   1824  CG  PHE A 314      51.324  33.730  71.250  1.00 29.69           C  
ATOM   1825  CD1 PHE A 314      50.950  32.976  70.141  1.00 30.03           C  
ATOM   1826  CD2 PHE A 314      50.989  33.276  72.517  1.00 29.58           C  
ATOM   1827  CE1 PHE A 314      50.261  31.773  70.297  1.00 30.30           C  
ATOM   1828  CE2 PHE A 314      50.312  32.076  72.678  1.00 29.97           C  
ATOM   1829  CZ  PHE A 314      49.942  31.326  71.565  1.00 29.66           C  
ATOM   1830  N   ARG A 315      50.049  36.860  72.707  1.00 29.46           N  
ATOM   1831  CA  ARG A 315      48.908  37.077  73.606  1.00 29.04           C  
ATOM   1832  C   ARG A 315      47.794  37.845  72.891  1.00 28.78           C  
ATOM   1833  O   ARG A 315      46.611  37.517  73.038  1.00 28.52           O  
ATOM   1834  CB  ARG A 315      49.332  37.820  74.873  1.00 28.82           C  
ATOM   1835  CG  ARG A 315      50.246  37.018  75.792  1.00 29.53           C  
ATOM   1836  CD  ARG A 315      50.178  37.570  77.199  1.00 29.63           C  
ATOM   1837  NE  ARG A 315      51.162  36.996  78.103  1.00 29.13           N  
ATOM   1838  CZ  ARG A 315      51.064  35.795  78.665  1.00 30.03           C  
ATOM   1839  NH1 ARG A 315      50.028  35.010  78.400  1.00 31.12           N  
ATOM   1840  NH2 ARG A 315      52.013  35.371  79.487  1.00 29.20           N  
ATOM   1841  N   GLN A 316      48.187  38.847  72.103  1.00 28.28           N  
ATOM   1842  CA  GLN A 316      47.234  39.674  71.360  1.00 28.42           C  
ATOM   1843  C   GLN A 316      46.494  38.857  70.298  1.00 28.06           C  
ATOM   1844  O   GLN A 316      45.290  39.045  70.092  1.00 27.65           O  
ATOM   1845  CB  GLN A 316      47.919  40.914  70.760  1.00 28.29           C  
ATOM   1846  CG  GLN A 316      46.997  41.877  70.017  1.00 29.95           C  
ATOM   1847  CD  GLN A 316      45.836  42.423  70.864  1.00 32.29           C  
ATOM   1848  OE1 GLN A 316      44.667  42.323  70.474  1.00 32.38           O  
ATOM   1849  NE2 GLN A 316      46.157  43.007  72.016  1.00 32.08           N  
ATOM   1850  N   PHE A 317      47.209  37.934  69.657  1.00 27.75           N  
ATOM   1851  CA  PHE A 317      46.585  37.009  68.713  1.00 27.81           C  
ATOM   1852  C   PHE A 317      45.476  36.196  69.408  1.00 27.61           C  
ATOM   1853  O   PHE A 317      44.369  36.049  68.879  1.00 27.54           O  
ATOM   1854  CB  PHE A 317      47.651  36.114  68.045  1.00 27.77           C  
ATOM   1855  CG  PHE A 317      47.098  34.860  67.404  1.00 27.96           C  
ATOM   1856  CD1 PHE A 317      46.019  34.913  66.529  1.00 27.97           C  
ATOM   1857  CD2 PHE A 317      47.676  33.624  67.670  1.00 28.69           C  
ATOM   1858  CE1 PHE A 317      45.515  33.752  65.948  1.00 28.75           C  
ATOM   1859  CE2 PHE A 317      47.179  32.453  67.090  1.00 28.06           C  
ATOM   1860  CZ  PHE A 317      46.104  32.517  66.228  1.00 28.09           C  
ATOM   1861  N   GLN A 318      45.778  35.697  70.604  1.00 27.47           N  
ATOM   1862  CA  GLN A 318      44.804  34.962  71.409  1.00 27.06           C  
ATOM   1863  C   GLN A 318      43.615  35.835  71.769  1.00 26.51           C  
ATOM   1864  O   GLN A 318      42.475  35.392  71.658  1.00 26.32           O  
ATOM   1865  CB  GLN A 318      45.450  34.417  72.680  1.00 27.17           C  
ATOM   1866  CG  GLN A 318      46.685  33.569  72.434  1.00 28.10           C  
ATOM   1867  CD  GLN A 318      47.126  32.851  73.675  1.00 28.93           C  
ATOM   1868  OE1 GLN A 318      46.523  31.858  74.064  1.00 30.67           O  
ATOM   1869  NE2 GLN A 318      48.181  33.346  74.311  1.00 28.37           N  
ATOM   1870  N   ARG A 319      43.884  37.072  72.189  1.00 25.95           N  
ATOM   1871  CA  ARG A 319      42.819  38.025  72.479  1.00 25.69           C  
ATOM   1872  C   ARG A 319      41.853  38.082  71.298  1.00 25.20           C  
ATOM   1873  O   ARG A 319      40.665  37.779  71.439  1.00 25.05           O  
ATOM   1874  CB  ARG A 319      43.379  39.427  72.729  1.00 25.87           C  
ATOM   1875  CG  ARG A 319      44.027  39.680  74.080  1.00 26.62           C  
ATOM   1876  CD  ARG A 319      44.362  41.170  74.200  1.00 28.00           C  
ATOM   1877  NE  ARG A 319      45.279  41.493  75.296  1.00 30.12           N  
ATOM   1878  CZ  ARG A 319      46.600  41.671  75.181  1.00 30.64           C  
ATOM   1879  NH1 ARG A 319      47.214  41.553  74.012  1.00 31.22           N  
ATOM   1880  NH2 ARG A 319      47.321  41.964  76.254  1.00 30.48           N  
ATOM   1881  N   ASN A 320      42.383  38.466  70.137  1.00 24.69           N  
ATOM   1882  CA  ASN A 320      41.594  38.655  68.928  1.00 24.40           C  
ATOM   1883  C   ASN A 320      40.823  37.396  68.551  1.00 24.79           C  
ATOM   1884  O   ASN A 320      39.660  37.474  68.145  1.00 24.80           O  
ATOM   1885  CB  ASN A 320      42.489  39.081  67.758  1.00 24.11           C  
ATOM   1886  CG  ASN A 320      43.161  40.423  67.982  1.00 23.44           C  
ATOM   1887  OD1 ASN A 320      42.753  41.204  68.843  1.00 24.97           O  
ATOM   1888  ND2 ASN A 320      44.194  40.700  67.201  1.00 20.57           N  
ATOM   1889  N   ALA A 321      41.475  36.238  68.688  1.00 24.79           N  
ATOM   1890  CA  ALA A 321      40.840  34.967  68.379  1.00 24.84           C  
ATOM   1891  C   ALA A 321      39.680  34.668  69.337  1.00 25.11           C  
ATOM   1892  O   ALA A 321      38.613  34.228  68.899  1.00 25.27           O  
ATOM   1893  CB  ALA A 321      41.861  33.843  68.386  1.00 24.52           C  
ATOM   1894  N   GLN A 322      39.899  34.912  70.632  1.00 25.15           N  
ATOM   1895  CA  GLN A 322      38.897  34.688  71.669  1.00 25.50           C  
ATOM   1896  C   GLN A 322      37.655  35.527  71.406  1.00 25.53           C  
ATOM   1897  O   GLN A 322      36.524  35.058  71.590  1.00 25.61           O  
ATOM   1898  CB  GLN A 322      39.470  35.038  73.044  1.00 25.57           C  
ATOM   1899  CG  GLN A 322      38.712  34.437  74.229  1.00 26.03           C  
ATOM   1900  CD  GLN A 322      39.435  34.636  75.557  1.00 26.40           C  
ATOM   1901  OE1 GLN A 322      39.082  34.022  76.561  1.00 27.56           O  
ATOM   1902  NE2 GLN A 322      40.444  35.505  75.570  1.00 27.44           N  
ATOM   1903  N   TYR A 323      37.872  36.762  70.968  1.00 25.36           N  
ATOM   1904  CA  TYR A 323      36.770  37.644  70.639  1.00 25.73           C  
ATOM   1905  C   TYR A 323      35.854  37.065  69.555  1.00 25.54           C  
ATOM   1906  O   TYR A 323      34.631  37.084  69.695  1.00 25.92           O  
ATOM   1907  CB  TYR A 323      37.269  39.022  70.200  1.00 26.32           C  
ATOM   1908  CG  TYR A 323      36.125  39.971  69.958  1.00 27.11           C  
ATOM   1909  CD1 TYR A 323      35.466  40.566  71.027  1.00 27.78           C  
ATOM   1910  CD2 TYR A 323      35.673  40.250  68.667  1.00 27.93           C  
ATOM   1911  CE1 TYR A 323      34.399  41.424  70.824  1.00 28.10           C  
ATOM   1912  CE2 TYR A 323      34.595  41.118  68.455  1.00 27.40           C  
ATOM   1913  CZ  TYR A 323      33.967  41.695  69.542  1.00 27.32           C  
ATOM   1914  OH  TYR A 323      32.904  42.561  69.372  1.00 28.22           O  
ATOM   1915  N   VAL A 324      36.446  36.561  68.481  1.00 24.78           N  
ATOM   1916  CA  VAL A 324      35.668  36.064  67.360  1.00 24.45           C  
ATOM   1917  C   VAL A 324      34.925  34.797  67.775  1.00 24.28           C  
ATOM   1918  O   VAL A 324      33.719  34.671  67.543  1.00 24.35           O  
ATOM   1919  CB  VAL A 324      36.559  35.836  66.111  1.00 24.26           C  
ATOM   1920  CG1 VAL A 324      35.776  35.157  64.991  1.00 23.96           C  
ATOM   1921  CG2 VAL A 324      37.134  37.167  65.639  1.00 23.68           C  
ATOM   1922  N   ALA A 325      35.648  33.891  68.429  1.00 23.93           N  
ATOM   1923  CA  ALA A 325      35.084  32.645  68.914  1.00 23.62           C  
ATOM   1924  C   ALA A 325      33.938  32.881  69.897  1.00 23.54           C  
ATOM   1925  O   ALA A 325      32.989  32.107  69.931  1.00 23.75           O  
ATOM   1926  CB  ALA A 325      36.169  31.775  69.544  1.00 23.44           C  
ATOM   1927  N   GLY A 326      34.014  33.956  70.675  1.00 23.33           N  
ATOM   1928  CA  GLY A 326      33.016  34.220  71.706  1.00 23.23           C  
ATOM   1929  C   GLY A 326      31.697  34.820  71.257  1.00 23.21           C  
ATOM   1930  O   GLY A 326      30.865  35.166  72.087  1.00 23.26           O  
ATOM   1931  N   LEU A 327      31.497  34.962  69.952  1.00 23.42           N  
ATOM   1932  CA  LEU A 327      30.254  35.548  69.440  1.00 23.48           C  
ATOM   1933  C   LEU A 327      29.211  34.461  69.171  1.00 23.40           C  
ATOM   1934  O   LEU A 327      28.020  34.657  69.405  1.00 23.20           O  
ATOM   1935  CB  LEU A 327      30.512  36.365  68.166  1.00 23.55           C  
ATOM   1936  CG  LEU A 327      31.565  37.486  68.127  1.00 23.64           C  
ATOM   1937  CD1 LEU A 327      31.795  37.956  66.693  1.00 23.74           C  
ATOM   1938  CD2 LEU A 327      31.194  38.666  68.997  1.00 22.84           C  
TER    1939      LEU A 327                                                      
ATOM   1940  N   HIS B  80      34.842  63.241  18.959  1.00 34.03           N  
ATOM   1941  CA  HIS B  80      34.810  61.884  18.333  1.00 33.90           C  
ATOM   1942  C   HIS B  80      35.464  60.809  19.203  1.00 34.05           C  
ATOM   1943  O   HIS B  80      35.025  59.654  19.185  1.00 34.23           O  
ATOM   1944  CB  HIS B  80      35.450  61.913  16.955  1.00 33.93           C  
ATOM   1945  N   ARG B  81      36.509  61.174  19.953  1.00 33.81           N  
ATOM   1946  CA  ARG B  81      37.147  60.230  20.883  1.00 33.46           C  
ATOM   1947  C   ARG B  81      36.168  59.866  21.991  1.00 33.25           C  
ATOM   1948  O   ARG B  81      36.173  58.735  22.479  1.00 33.64           O  
ATOM   1949  CB  ARG B  81      38.436  60.802  21.469  1.00 33.64           C  
ATOM   1950  N   VAL B  82      35.331  60.833  22.372  1.00 32.59           N  
ATOM   1951  CA  VAL B  82      34.206  60.618  23.292  1.00 31.85           C  
ATOM   1952  C   VAL B  82      33.324  59.447  22.838  1.00 31.01           C  
ATOM   1953  O   VAL B  82      33.039  58.534  23.617  1.00 31.20           O  
ATOM   1954  CB  VAL B  82      33.308  61.898  23.408  1.00 32.21           C  
ATOM   1955  CG1 VAL B  82      32.164  61.682  24.418  1.00 32.28           C  
ATOM   1956  CG2 VAL B  82      34.144  63.140  23.764  1.00 31.98           C  
ATOM   1957  N   THR B  83      32.909  59.491  21.571  1.00 29.94           N  
ATOM   1958  CA  THR B  83      31.964  58.533  20.990  1.00 28.85           C  
ATOM   1959  C   THR B  83      32.597  57.155  20.777  1.00 28.03           C  
ATOM   1960  O   THR B  83      31.966  56.128  21.034  1.00 27.89           O  
ATOM   1961  CB  THR B  83      31.387  59.061  19.653  1.00 28.86           C  
ATOM   1962  OG1 THR B  83      31.039  60.450  19.789  1.00 29.03           O  
ATOM   1963  CG2 THR B  83      30.152  58.262  19.237  1.00 28.54           C  
ATOM   1964  N   LEU B  84      33.842  57.143  20.312  1.00 27.07           N  
ATOM   1965  CA  LEU B  84      34.597  55.910  20.164  1.00 26.17           C  
ATOM   1966  C   LEU B  84      34.680  55.188  21.506  1.00 25.86           C  
ATOM   1967  O   LEU B  84      34.384  54.002  21.599  1.00 25.78           O  
ATOM   1968  CB  LEU B  84      35.997  56.200  19.620  1.00 26.08           C  
ATOM   1969  CG  LEU B  84      36.604  55.228  18.599  1.00 25.70           C  
ATOM   1970  CD1 LEU B  84      38.104  55.456  18.436  1.00 23.69           C  
ATOM   1971  CD2 LEU B  84      36.318  53.772  18.956  1.00 24.92           C  
ATOM   1972  N   ARG B  85      35.066  55.924  22.546  1.00 25.54           N  
ATOM   1973  CA  ARG B  85      35.164  55.392  23.907  1.00 24.94           C  
ATOM   1974  C   ARG B  85      33.847  54.747  24.385  1.00 24.65           C  
ATOM   1975  O   ARG B  85      33.836  53.584  24.782  1.00 24.36           O  
ATOM   1976  CB  ARG B  85      35.614  56.510  24.850  1.00 24.88           C  
ATOM   1977  CG  ARG B  85      35.967  56.064  26.251  1.00 24.45           C  
ATOM   1978  CD  ARG B  85      36.102  57.265  27.160  1.00 23.68           C  
ATOM   1979  NE  ARG B  85      36.122  56.868  28.558  1.00 24.12           N  
ATOM   1980  CZ  ARG B  85      35.038  56.701  29.310  1.00 23.97           C  
ATOM   1981  NH1 ARG B  85      33.824  56.901  28.810  1.00 23.83           N  
ATOM   1982  NH2 ARG B  85      35.173  56.322  30.569  1.00 24.46           N  
ATOM   1983  N   LYS B  86      32.748  55.501  24.327  1.00 24.43           N  
ATOM   1984  CA  LYS B  86      31.420  54.996  24.677  1.00 24.65           C  
ATOM   1985  C   LYS B  86      31.131  53.678  23.947  1.00 24.46           C  
ATOM   1986  O   LYS B  86      30.735  52.690  24.566  1.00 24.64           O  
ATOM   1987  CB  LYS B  86      30.353  56.049  24.344  1.00 24.62           C  
ATOM   1988  CG  LYS B  86      28.932  55.741  24.818  1.00 24.69           C  
ATOM   1989  CD  LYS B  86      28.047  56.977  24.642  1.00 26.08           C  
ATOM   1990  CE  LYS B  86      26.572  56.704  24.956  1.00 29.27           C  
ATOM   1991  NZ  LYS B  86      26.319  56.498  26.427  1.00 31.01           N  
ATOM   1992  N   ALA B  87      31.364  53.675  22.637  1.00 24.16           N  
ATOM   1993  CA  ALA B  87      31.176  52.503  21.785  1.00 23.76           C  
ATOM   1994  C   ALA B  87      32.060  51.307  22.161  1.00 23.36           C  
ATOM   1995  O   ALA B  87      31.600  50.172  22.111  1.00 23.22           O  
ATOM   1996  CB  ALA B  87      31.385  52.888  20.320  1.00 23.66           C  
ATOM   1997  N   THR B  88      33.313  51.563  22.536  1.00 23.32           N  
ATOM   1998  CA  THR B  88      34.249  50.501  22.921  1.00 23.51           C  
ATOM   1999  C   THR B  88      33.813  49.865  24.237  1.00 23.95           C  
ATOM   2000  O   THR B  88      33.753  48.644  24.348  1.00 24.22           O  
ATOM   2001  CB  THR B  88      35.702  51.021  23.040  1.00 23.35           C  
ATOM   2002  OG1 THR B  88      36.157  51.472  21.757  1.00 24.21           O  
ATOM   2003  CG2 THR B  88      36.647  49.932  23.546  1.00 22.24           C  
ATOM   2004  N   LEU B  89      33.501  50.696  25.228  1.00 24.30           N  
ATOM   2005  CA  LEU B  89      33.018  50.199  26.517  1.00 24.40           C  
ATOM   2006  C   LEU B  89      31.669  49.490  26.365  1.00 24.43           C  
ATOM   2007  O   LEU B  89      31.428  48.470  27.015  1.00 24.52           O  
ATOM   2008  CB  LEU B  89      32.913  51.333  27.553  1.00 24.29           C  
ATOM   2009  CG  LEU B  89      34.163  52.085  28.038  1.00 23.72           C  
ATOM   2010  CD1 LEU B  89      33.751  53.325  28.806  1.00 22.22           C  
ATOM   2011  CD2 LEU B  89      35.054  51.204  28.899  1.00 23.11           C  
ATOM   2012  N   ALA B  90      30.802  50.020  25.504  1.00 24.37           N  
ATOM   2013  CA  ALA B  90      29.479  49.431  25.296  1.00 24.80           C  
ATOM   2014  C   ALA B  90      29.552  48.052  24.653  1.00 25.20           C  
ATOM   2015  O   ALA B  90      28.841  47.136  25.079  1.00 25.18           O  
ATOM   2016  CB  ALA B  90      28.579  50.355  24.480  1.00 24.76           C  
ATOM   2017  N   SER B  91      30.414  47.898  23.647  1.00 25.58           N  
ATOM   2018  CA  SER B  91      30.526  46.615  22.948  1.00 26.16           C  
ATOM   2019  C   SER B  91      31.206  45.552  23.815  1.00 26.58           C  
ATOM   2020  O   SER B  91      30.885  44.359  23.712  1.00 26.65           O  
ATOM   2021  CB  SER B  91      31.235  46.766  21.600  1.00 26.14           C  
ATOM   2022  OG  SER B  91      32.434  46.019  21.539  1.00 26.37           O  
ATOM   2023  N   LEU B  92      32.137  45.992  24.660  1.00 26.79           N  
ATOM   2024  CA  LEU B  92      32.769  45.127  25.645  1.00 27.44           C  
ATOM   2025  C   LEU B  92      31.765  44.610  26.674  1.00 28.25           C  
ATOM   2026  O   LEU B  92      31.683  43.410  26.905  1.00 28.33           O  
ATOM   2027  CB  LEU B  92      33.911  45.858  26.341  1.00 27.11           C  
ATOM   2028  CG  LEU B  92      34.584  45.156  27.515  1.00 26.20           C  
ATOM   2029  CD1 LEU B  92      35.327  43.900  27.063  1.00 25.08           C  
ATOM   2030  CD2 LEU B  92      35.524  46.146  28.175  1.00 24.86           C  
HETATM 2031  N   MSE B  93      30.996  45.513  27.271  1.00 29.59           N  
HETATM 2032  CA  MSE B  93      29.937  45.131  28.208  1.00 31.32           C  
HETATM 2033  C   MSE B  93      28.944  44.142  27.600  1.00 31.44           C  
HETATM 2034  O   MSE B  93      28.546  43.179  28.264  1.00 31.63           O  
HETATM 2035  CB  MSE B  93      29.177  46.358  28.726  1.00 32.38           C  
HETATM 2036  CG  MSE B  93      29.958  47.267  29.681  1.00 34.83           C  
HETATM 2037 SE   MSE B  93      30.788  46.350  31.205  1.00 42.85          SE  
HETATM 2038  CE  MSE B  93      32.654  46.284  30.602  1.00 36.60           C  
ATOM   2039  N   GLN B  94      28.552  44.383  26.347  1.00 31.19           N  
ATOM   2040  CA  GLN B  94      27.632  43.502  25.628  1.00 31.39           C  
ATOM   2041  C   GLN B  94      28.158  42.072  25.495  1.00 30.78           C  
ATOM   2042  O   GLN B  94      27.393  41.112  25.628  1.00 30.95           O  
ATOM   2043  CB  GLN B  94      27.311  44.069  24.241  1.00 31.47           C  
ATOM   2044  CG  GLN B  94      26.003  44.843  24.167  1.00 32.69           C  
ATOM   2045  CD  GLN B  94      25.823  45.592  22.849  1.00 32.90           C  
ATOM   2046  OE1 GLN B  94      24.896  45.308  22.073  1.00 33.99           O  
ATOM   2047  NE2 GLN B  94      26.711  46.556  22.589  1.00 34.83           N  
ATOM   2048  N   SER B  95      29.457  41.934  25.243  1.00 30.19           N  
ATOM   2049  CA  SER B  95      30.075  40.618  25.054  1.00 29.82           C  
ATOM   2050  C   SER B  95      30.125  39.768  26.331  1.00 29.75           C  
ATOM   2051  O   SER B  95      30.374  38.568  26.274  1.00 29.84           O  
ATOM   2052  CB  SER B  95      31.488  40.768  24.493  1.00 29.73           C  
ATOM   2053  OG  SER B  95      32.407  41.170  25.500  1.00 29.20           O  
ATOM   2054  N   LEU B  96      29.880  40.395  27.476  1.00 29.60           N  
ATOM   2055  CA  LEU B  96      30.081  39.746  28.766  1.00 29.38           C  
ATOM   2056  C   LEU B  96      28.953  38.806  29.163  1.00 29.32           C  
ATOM   2057  O   LEU B  96      29.039  38.133  30.193  1.00 29.63           O  
ATOM   2058  CB  LEU B  96      30.341  40.793  29.852  1.00 29.50           C  
ATOM   2059  CG  LEU B  96      31.787  41.211  30.165  1.00 29.28           C  
ATOM   2060  CD1 LEU B  96      32.767  41.067  29.008  1.00 28.18           C  
ATOM   2061  CD2 LEU B  96      31.799  42.623  30.704  1.00 29.25           C  
ATOM   2062  N   SER B  97      27.905  38.752  28.345  1.00 29.14           N  
ATOM   2063  CA  SER B  97      26.869  37.729  28.478  1.00 28.86           C  
ATOM   2064  C   SER B  97      27.330  36.446  27.765  1.00 28.73           C  
ATOM   2065  O   SER B  97      26.681  35.390  27.850  1.00 28.46           O  
ATOM   2066  CB  SER B  97      25.554  38.236  27.897  1.00 28.85           C  
ATOM   2067  OG  SER B  97      25.692  38.440  26.504  1.00 29.52           O  
ATOM   2068  N   GLY B  98      28.451  36.566  27.052  1.00 28.61           N  
ATOM   2069  CA  GLY B  98      29.167  35.427  26.481  1.00 28.44           C  
ATOM   2070  C   GLY B  98      30.121  34.818  27.490  1.00 28.44           C  
ATOM   2071  O   GLY B  98      30.904  35.524  28.143  1.00 28.43           O  
ATOM   2072  N   GLU B  99      30.037  33.495  27.604  1.00 28.54           N  
ATOM   2073  CA  GLU B  99      30.787  32.667  28.561  1.00 28.52           C  
ATOM   2074  C   GLU B  99      32.290  32.986  28.627  1.00 28.32           C  
ATOM   2075  O   GLU B  99      32.805  33.402  29.675  1.00 28.31           O  
ATOM   2076  CB  GLU B  99      30.560  31.181  28.216  1.00 28.30           C  
ATOM   2077  CG  GLU B  99      30.913  30.190  29.316  1.00 28.74           C  
ATOM   2078  CD  GLU B  99      30.846  28.732  28.863  1.00 28.84           C  
ATOM   2079  OE1 GLU B  99      29.988  28.396  28.021  1.00 30.03           O  
ATOM   2080  OE2 GLU B  99      31.646  27.909  29.357  1.00 29.43           O  
ATOM   2081  N   SER B 100      32.974  32.788  27.501  1.00 28.13           N  
ATOM   2082  CA  SER B 100      34.423  32.955  27.394  1.00 28.09           C  
ATOM   2083  C   SER B 100      34.859  34.393  27.717  1.00 28.04           C  
ATOM   2084  O   SER B 100      35.822  34.610  28.458  1.00 28.10           O  
ATOM   2085  CB  SER B 100      34.901  32.549  25.987  1.00 28.37           C  
ATOM   2086  OG  SER B 100      34.309  31.330  25.552  1.00 28.00           O  
ATOM   2087  N   SER B 101      34.147  35.367  27.154  1.00 27.78           N  
ATOM   2088  CA  SER B 101      34.426  36.770  27.409  1.00 27.49           C  
ATOM   2089  C   SER B 101      34.292  37.052  28.896  1.00 27.45           C  
ATOM   2090  O   SER B 101      35.213  37.587  29.513  1.00 27.33           O  
ATOM   2091  CB  SER B 101      33.474  37.666  26.617  1.00 27.64           C  
ATOM   2092  OG  SER B 101      33.386  37.251  25.270  1.00 28.02           O  
ATOM   2093  N   ASN B 102      33.153  36.660  29.466  1.00 27.52           N  
ATOM   2094  CA  ASN B 102      32.876  36.871  30.885  1.00 27.67           C  
ATOM   2095  C   ASN B 102      33.993  36.377  31.794  1.00 27.52           C  
ATOM   2096  O   ASN B 102      34.365  37.060  32.752  1.00 27.59           O  
ATOM   2097  CB  ASN B 102      31.543  36.232  31.283  1.00 27.78           C  
ATOM   2098  CG  ASN B 102      31.061  36.694  32.643  1.00 28.93           C  
ATOM   2099  OD1 ASN B 102      30.187  37.552  32.742  1.00 30.93           O  
ATOM   2100  ND2 ASN B 102      31.642  36.144  33.701  1.00 29.71           N  
ATOM   2101  N   ARG B 103      34.540  35.204  31.481  1.00 27.61           N  
ATOM   2102  CA  ARG B 103      35.530  34.561  32.352  1.00 27.58           C  
ATOM   2103  C   ARG B 103      36.931  35.124  32.169  1.00 27.73           C  
ATOM   2104  O   ARG B 103      37.760  35.062  33.072  1.00 27.44           O  
ATOM   2105  CB  ARG B 103      35.544  33.051  32.123  1.00 27.73           C  
ATOM   2106  CG  ARG B 103      34.307  32.343  32.653  1.00 26.86           C  
ATOM   2107  CD  ARG B 103      34.570  30.874  32.900  1.00 24.91           C  
ATOM   2108  NE  ARG B 103      33.317  30.148  33.066  1.00 23.50           N  
ATOM   2109  CZ  ARG B 103      32.814  29.306  32.172  1.00 22.52           C  
ATOM   2110  NH1 ARG B 103      33.466  29.055  31.044  1.00 22.13           N  
ATOM   2111  NH2 ARG B 103      31.661  28.703  32.416  1.00 22.77           N  
ATOM   2112  N   VAL B 104      37.192  35.654  30.982  1.00 27.99           N  
ATOM   2113  CA  VAL B 104      38.462  36.282  30.699  1.00 28.26           C  
ATOM   2114  C   VAL B 104      38.485  37.628  31.391  1.00 28.83           C  
ATOM   2115  O   VAL B 104      39.515  38.029  31.944  1.00 29.11           O  
ATOM   2116  CB  VAL B 104      38.692  36.434  29.178  1.00 28.32           C  
ATOM   2117  CG1 VAL B 104      39.814  37.431  28.885  1.00 27.49           C  
ATOM   2118  CG2 VAL B 104      39.005  35.076  28.565  1.00 27.67           C  
HETATM 2119  N   MSE B 105      37.341  38.284  31.397  1.00 20.00           N  
HETATM 2120  CA  MSE B 105      37.180  39.596  31.961  1.00 20.00           C  
HETATM 2121  C   MSE B 105      37.200  39.653  33.470  1.00 20.00           C  
HETATM 2122  O   MSE B 105      37.846  40.476  34.022  1.00 29.65           O  
HETATM 2123  CB  MSE B 105      35.922  40.250  31.412  1.00 20.00           C  
HETATM 2124  CG  MSE B 105      35.748  41.712  31.729  1.00 20.00           C  
HETATM 2125 SE   MSE B 105      36.882  42.863  30.836  1.00 20.00          SE  
HETATM 2126  CE  MSE B 105      37.969  43.425  32.189  1.00 20.00           C  
ATOM   2127  N   TRP B 106      36.488  38.769  34.131  1.00 29.45           N  
ATOM   2128  CA  TRP B 106      36.372  38.831  35.599  1.00 29.67           C  
ATOM   2129  C   TRP B 106      37.176  37.786  36.382  1.00 29.90           C  
ATOM   2130  O   TRP B 106      36.808  37.406  37.483  1.00 30.42           O  
ATOM   2131  CB  TRP B 106      34.896  38.845  36.041  1.00 29.61           C  
ATOM   2132  CG  TRP B 106      34.007  39.850  35.293  1.00 29.82           C  
ATOM   2133  CD1 TRP B 106      32.822  39.580  34.672  1.00 29.91           C  
ATOM   2134  CD2 TRP B 106      34.246  41.257  35.096  1.00 29.49           C  
ATOM   2135  NE1 TRP B 106      32.308  40.720  34.102  1.00 30.20           N  
ATOM   2136  CE2 TRP B 106      33.157  41.765  34.347  1.00 30.09           C  
ATOM   2137  CE3 TRP B 106      35.263  42.137  35.483  1.00 29.45           C  
ATOM   2138  CZ2 TRP B 106      33.056  43.116  33.980  1.00 28.74           C  
ATOM   2139  CZ3 TRP B 106      35.165  43.482  35.107  1.00 29.91           C  
ATOM   2140  CH2 TRP B 106      34.067  43.953  34.363  1.00 28.71           C  
ATOM   2141  N   ASN B 107      38.278  37.323  35.814  1.00 30.25           N  
ATOM   2142  CA  ASN B 107      39.211  36.477  36.528  1.00 30.43           C  
ATOM   2143  C   ASN B 107      39.986  37.346  37.502  1.00 30.79           C  
ATOM   2144  O   ASN B 107      40.115  38.549  37.281  1.00 30.45           O  
ATOM   2145  CB  ASN B 107      40.177  35.820  35.544  1.00 30.69           C  
ATOM   2146  CG  ASN B 107      40.850  34.594  36.121  1.00 30.90           C  
ATOM   2147  OD1 ASN B 107      42.050  34.598  36.414  1.00 31.09           O  
ATOM   2148  ND2 ASN B 107      40.074  33.533  36.296  1.00 31.28           N  
ATOM   2149  N   ASP B 108      40.499  36.747  38.577  1.00 31.54           N  
ATOM   2150  CA  ASP B 108      41.218  37.515  39.594  1.00 32.17           C  
ATOM   2151  C   ASP B 108      42.718  37.217  39.679  1.00 32.69           C  
ATOM   2152  O   ASP B 108      43.358  37.512  40.695  1.00 32.80           O  
ATOM   2153  CB  ASP B 108      40.541  37.404  40.969  1.00 32.40           C  
ATOM   2154  CG  ASP B 108      40.497  35.977  41.503  1.00 32.79           C  
ATOM   2155  OD1 ASP B 108      40.783  35.034  40.735  1.00 33.95           O  
ATOM   2156  OD2 ASP B 108      40.159  35.805  42.697  1.00 31.74           O  
ATOM   2157  N   ARG B 109      43.273  36.636  38.617  1.00 33.11           N  
ATOM   2158  CA  ARG B 109      44.724  36.553  38.480  1.00 33.69           C  
ATOM   2159  C   ARG B 109      45.210  37.895  37.942  1.00 34.04           C  
ATOM   2160  O   ARG B 109      44.721  38.380  36.916  1.00 34.18           O  
ATOM   2161  CB  ARG B 109      45.140  35.395  37.559  1.00 33.76           C  
ATOM   2162  N   TYR B 110      46.154  38.508  38.646  1.00 34.34           N  
ATOM   2163  CA  TYR B 110      46.586  39.859  38.304  1.00 34.41           C  
ATOM   2164  C   TYR B 110      48.048  39.925  37.867  1.00 34.47           C  
ATOM   2165  O   TYR B 110      48.430  40.831  37.121  1.00 34.55           O  
ATOM   2166  CB  TYR B 110      46.343  40.812  39.480  1.00 34.45           C  
ATOM   2167  CG  TYR B 110      44.908  40.900  39.957  1.00 34.37           C  
ATOM   2168  CD1 TYR B 110      43.860  41.119  39.057  1.00 34.35           C  
ATOM   2169  CD2 TYR B 110      44.600  40.795  41.312  1.00 34.21           C  
ATOM   2170  CE1 TYR B 110      42.538  41.203  39.493  1.00 34.29           C  
ATOM   2171  CE2 TYR B 110      43.283  40.890  41.764  1.00 34.48           C  
ATOM   2172  CZ  TYR B 110      42.254  41.095  40.849  1.00 34.95           C  
ATOM   2173  OH  TYR B 110      40.945  41.192  41.289  1.00 35.09           O  
ATOM   2174  N   ASP B 111      48.858  38.968  38.326  1.00 34.34           N  
ATOM   2175  CA  ASP B 111      50.292  38.948  38.006  1.00 34.28           C  
ATOM   2176  C   ASP B 111      50.591  38.449  36.581  1.00 34.13           C  
ATOM   2177  O   ASP B 111      51.752  38.258  36.198  1.00 34.38           O  
ATOM   2178  CB  ASP B 111      51.109  38.197  39.090  1.00 34.37           C  
ATOM   2179  CG  ASP B 111      50.913  36.676  39.067  1.00 34.17           C  
ATOM   2180  OD1 ASP B 111      49.783  36.198  38.819  1.00 33.41           O  
ATOM   2181  OD2 ASP B 111      51.910  35.958  39.328  1.00 34.10           O  
ATOM   2182  N   THR B 112      49.525  38.275  35.803  1.00 33.73           N  
ATOM   2183  CA  THR B 112      49.592  37.786  34.433  1.00 33.34           C  
ATOM   2184  C   THR B 112      48.860  38.760  33.517  1.00 32.67           C  
ATOM   2185  O   THR B 112      47.835  39.333  33.893  1.00 32.48           O  
ATOM   2186  CB  THR B 112      48.940  36.382  34.324  1.00 33.64           C  
ATOM   2187  OG1 THR B 112      49.533  35.505  35.292  1.00 34.45           O  
ATOM   2188  CG2 THR B 112      49.109  35.782  32.923  1.00 33.57           C  
ATOM   2189  N   LEU B 113      49.399  38.941  32.317  1.00 31.99           N  
ATOM   2190  CA  LEU B 113      48.801  39.798  31.303  1.00 31.40           C  
ATOM   2191  C   LEU B 113      47.364  39.364  31.035  1.00 30.88           C  
ATOM   2192  O   LEU B 113      47.108  38.195  30.762  1.00 30.86           O  
ATOM   2193  CB  LEU B 113      49.640  39.707  30.031  1.00 31.44           C  
ATOM   2194  CG  LEU B 113      49.574  40.793  28.966  1.00 31.46           C  
ATOM   2195  CD1 LEU B 113      49.933  42.161  29.541  1.00 31.36           C  
ATOM   2196  CD2 LEU B 113      50.514  40.416  27.829  1.00 31.30           C  
ATOM   2197  N   LEU B 114      46.428  40.303  31.138  1.00 30.44           N  
ATOM   2198  CA  LEU B 114      45.003  39.984  31.021  1.00 30.06           C  
ATOM   2199  C   LEU B 114      44.632  39.285  29.718  1.00 29.94           C  
ATOM   2200  O   LEU B 114      43.796  38.389  29.723  1.00 30.19           O  
ATOM   2201  CB  LEU B 114      44.139  41.233  31.185  1.00 29.99           C  
ATOM   2202  CG  LEU B 114      42.638  40.993  31.422  1.00 30.00           C  
ATOM   2203  CD1 LEU B 114      42.055  42.085  32.341  1.00 30.23           C  
ATOM   2204  CD2 LEU B 114      41.835  40.874  30.130  1.00 28.45           C  
ATOM   2205  N   ILE B 115      45.236  39.699  28.608  1.00 29.63           N  
ATOM   2206  CA  ILE B 115      44.911  39.115  27.300  1.00 29.23           C  
ATOM   2207  C   ILE B 115      45.544  37.739  27.079  1.00 28.95           C  
ATOM   2208  O   ILE B 115      45.352  37.130  26.027  1.00 28.86           O  
ATOM   2209  CB  ILE B 115      45.261  40.070  26.125  1.00 29.34           C  
ATOM   2210  CG1 ILE B 115      46.784  40.297  26.040  1.00 29.33           C  
ATOM   2211  CG2 ILE B 115      44.451  41.373  26.251  1.00 29.19           C  
ATOM   2212  CD1 ILE B 115      47.240  41.288  24.963  1.00 29.33           C  
ATOM   2213  N   ALA B 116      46.298  37.261  28.071  1.00 28.62           N  
ATOM   2214  CA  ALA B 116      46.900  35.929  28.025  1.00 28.31           C  
ATOM   2215  C   ALA B 116      46.057  34.915  28.807  1.00 28.19           C  
ATOM   2216  O   ALA B 116      46.382  33.728  28.855  1.00 28.19           O  
ATOM   2217  CB  ALA B 116      48.331  35.970  28.548  1.00 28.12           C  
ATOM   2218  N   ARG B 117      44.973  35.395  29.410  1.00 27.97           N  
ATOM   2219  CA  ARG B 117      44.001  34.536  30.078  1.00 27.98           C  
ATOM   2220  C   ARG B 117      43.313  33.638  29.060  1.00 28.34           C  
ATOM   2221  O   ARG B 117      42.672  34.113  28.120  1.00 28.35           O  
ATOM   2222  CB  ARG B 117      42.965  35.371  30.835  1.00 27.71           C  
ATOM   2223  CG  ARG B 117      43.558  36.196  31.957  1.00 27.06           C  
ATOM   2224  CD  ARG B 117      42.514  37.071  32.595  1.00 26.25           C  
ATOM   2225  NE  ARG B 117      43.038  37.781  33.762  1.00 25.97           N  
ATOM   2226  CZ  ARG B 117      42.379  38.724  34.434  1.00 25.74           C  
ATOM   2227  NH1 ARG B 117      41.149  39.095  34.073  1.00 25.37           N  
ATOM   2228  NH2 ARG B 117      42.960  39.310  35.470  1.00 24.80           N  
ATOM   2229  N   ASP B 118      43.471  32.333  29.243  1.00 28.71           N  
ATOM   2230  CA  ASP B 118      42.868  31.369  28.351  1.00 29.02           C  
ATOM   2231  C   ASP B 118      41.522  30.902  28.915  1.00 29.21           C  
ATOM   2232  O   ASP B 118      41.471  30.368  30.028  1.00 29.25           O  
ATOM   2233  CB  ASP B 118      43.829  30.201  28.146  1.00 29.17           C  
ATOM   2234  CG  ASP B 118      43.265  29.113  27.257  1.00 29.52           C  
ATOM   2235  OD1 ASP B 118      42.427  29.395  26.376  1.00 29.44           O  
ATOM   2236  OD2 ASP B 118      43.685  27.957  27.440  1.00 31.25           O  
ATOM   2237  N   PRO B 119      40.428  31.116  28.150  1.00 29.45           N  
ATOM   2238  CA  PRO B 119      39.063  30.704  28.527  1.00 29.58           C  
ATOM   2239  C   PRO B 119      38.971  29.205  28.798  1.00 29.56           C  
ATOM   2240  O   PRO B 119      38.290  28.784  29.735  1.00 29.46           O  
ATOM   2241  CB  PRO B 119      38.229  31.045  27.283  1.00 29.57           C  
ATOM   2242  CG  PRO B 119      39.015  32.048  26.541  1.00 29.50           C  
ATOM   2243  CD  PRO B 119      40.456  31.792  26.839  1.00 29.33           C  
ATOM   2244  N   ARG B 120      39.648  28.436  27.947  1.00 29.58           N  
ATOM   2245  CA  ARG B 120      39.806  26.995  28.046  1.00 29.72           C  
ATOM   2246  C   ARG B 120      40.362  26.604  29.420  1.00 29.56           C  
ATOM   2247  O   ARG B 120      39.723  25.870  30.179  1.00 29.25           O  
ATOM   2248  CB  ARG B 120      40.805  26.586  26.974  1.00 30.05           C  
ATOM   2249  CG  ARG B 120      40.500  25.359  26.157  1.00 31.88           C  
ATOM   2250  CD  ARG B 120      41.682  25.137  25.208  1.00 34.01           C  
ATOM   2251  NE  ARG B 120      41.365  24.298  24.052  1.00 35.79           N  
ATOM   2252  CZ  ARG B 120      41.021  24.764  22.851  1.00 36.49           C  
ATOM   2253  NH1 ARG B 120      40.929  26.073  22.630  1.00 35.91           N  
ATOM   2254  NH2 ARG B 120      40.764  23.913  21.863  1.00 36.75           N  
ATOM   2255  N   GLU B 121      41.552  27.121  29.733  1.00 29.43           N  
ATOM   2256  CA  GLU B 121      42.269  26.825  30.984  1.00 29.31           C  
ATOM   2257  C   GLU B 121      41.461  27.206  32.234  1.00 28.56           C  
ATOM   2258  O   GLU B 121      41.517  26.515  33.254  1.00 28.09           O  
ATOM   2259  CB  GLU B 121      43.633  27.538  30.976  1.00 29.60           C  
ATOM   2260  CG  GLU B 121      44.615  27.120  32.077  1.00 31.92           C  
ATOM   2261  CD  GLU B 121      44.344  27.787  33.430  1.00 34.57           C  
ATOM   2262  OE1 GLU B 121      43.853  28.939  33.448  1.00 34.95           O  
ATOM   2263  OE2 GLU B 121      44.625  27.153  34.478  1.00 35.54           O  
ATOM   2264  N   ILE B 122      40.721  28.309  32.136  1.00 27.98           N  
ATOM   2265  CA  ILE B 122      39.899  28.806  33.232  1.00 27.28           C  
ATOM   2266  C   ILE B 122      38.709  27.881  33.509  1.00 27.16           C  
ATOM   2267  O   ILE B 122      38.363  27.651  34.670  1.00 27.20           O  
ATOM   2268  CB  ILE B 122      39.431  30.275  32.981  1.00 27.46           C  
ATOM   2269  CG1 ILE B 122      40.616  31.249  33.088  1.00 27.15           C  
ATOM   2270  CG2 ILE B 122      38.297  30.683  33.940  1.00 26.76           C  
ATOM   2271  CD1 ILE B 122      40.335  32.671  32.541  1.00 26.72           C  
ATOM   2272  N   LYS B 123      38.091  27.352  32.454  1.00 26.82           N  
ATOM   2273  CA  LYS B 123      36.953  26.447  32.618  1.00 26.75           C  
ATOM   2274  C   LYS B 123      37.397  25.107  33.207  1.00 26.22           C  
ATOM   2275  O   LYS B 123      36.698  24.521  34.040  1.00 26.15           O  
ATOM   2276  CB  LYS B 123      36.189  26.237  31.301  1.00 27.01           C  
ATOM   2277  CG  LYS B 123      34.885  25.472  31.498  1.00 28.18           C  
ATOM   2278  CD  LYS B 123      33.999  25.485  30.277  1.00 31.09           C  
ATOM   2279  CE  LYS B 123      32.667  24.786  30.573  1.00 32.02           C  
ATOM   2280  NZ  LYS B 123      31.819  24.695  29.340  1.00 33.09           N  
ATOM   2281  N   ASN B 124      38.560  24.633  32.779  1.00 25.85           N  
ATOM   2282  CA  ASN B 124      39.145  23.427  33.351  1.00 25.69           C  
ATOM   2283  C   ASN B 124      39.461  23.601  34.823  1.00 25.46           C  
ATOM   2284  O   ASN B 124      39.287  22.665  35.608  1.00 25.39           O  
ATOM   2285  CB  ASN B 124      40.386  22.995  32.580  1.00 25.74           C  
ATOM   2286  CG  ASN B 124      40.044  22.415  31.222  1.00 26.73           C  
ATOM   2287  OD1 ASN B 124      38.955  21.875  31.027  1.00 27.83           O  
ATOM   2288  ND2 ASN B 124      40.972  22.523  30.274  1.00 27.05           N  
ATOM   2289  N   ALA B 125      39.904  24.802  35.196  1.00 25.03           N  
ATOM   2290  CA  ALA B 125      40.188  25.114  36.592  1.00 25.04           C  
ATOM   2291  C   ALA B 125      38.910  24.994  37.441  1.00 25.38           C  
ATOM   2292  O   ALA B 125      38.932  24.386  38.508  1.00 25.20           O  
ATOM   2293  CB  ALA B 125      40.837  26.495  36.722  1.00 24.32           C  
ATOM   2294  N   ILE B 126      37.802  25.545  36.938  1.00 25.91           N  
ATOM   2295  CA  ILE B 126      36.492  25.438  37.584  1.00 26.63           C  
ATOM   2296  C   ILE B 126      36.038  23.984  37.679  1.00 27.10           C  
ATOM   2297  O   ILE B 126      35.654  23.535  38.754  1.00 26.86           O  
ATOM   2298  CB  ILE B 126      35.394  26.286  36.855  1.00 26.88           C  
ATOM   2299  CG1 ILE B 126      35.750  27.775  36.865  1.00 26.10           C  
ATOM   2300  CG2 ILE B 126      34.022  26.098  37.517  1.00 26.71           C  
ATOM   2301  CD1 ILE B 126      35.054  28.562  35.791  1.00 25.16           C  
ATOM   2302  N   GLU B 127      36.084  23.263  36.554  1.00 27.93           N  
ATOM   2303  CA  GLU B 127      35.767  21.830  36.525  1.00 28.85           C  
ATOM   2304  C   GLU B 127      36.550  21.032  37.581  1.00 28.78           C  
ATOM   2305  O   GLU B 127      35.958  20.245  38.318  1.00 29.03           O  
ATOM   2306  CB  GLU B 127      36.006  21.237  35.135  1.00 28.76           C  
ATOM   2307  CG  GLU B 127      34.897  21.504  34.105  1.00 29.79           C  
ATOM   2308  CD  GLU B 127      35.333  21.188  32.666  1.00 30.28           C  
ATOM   2309  OE1 GLU B 127      35.850  20.076  32.409  1.00 33.18           O  
ATOM   2310  OE2 GLU B 127      35.164  22.055  31.781  1.00 32.93           O  
ATOM   2311  N   LYS B 128      37.863  21.256  37.661  1.00 28.76           N  
ATOM   2312  CA  LYS B 128      38.728  20.567  38.626  1.00 28.88           C  
ATOM   2313  C   LYS B 128      38.383  20.935  40.069  1.00 29.38           C  
ATOM   2314  O   LYS B 128      38.314  20.058  40.944  1.00 29.63           O  
ATOM   2315  CB  LYS B 128      40.205  20.844  38.321  1.00 28.64           C  
ATOM   2316  CG  LYS B 128      41.157  20.494  39.431  1.00 28.45           C  
ATOM   2317  CD  LYS B 128      42.526  20.132  38.906  1.00 29.51           C  
ATOM   2318  CE  LYS B 128      43.455  19.686  40.037  1.00 30.44           C  
ATOM   2319  NZ  LYS B 128      42.715  18.940  41.113  1.00 30.32           N  
ATOM   2320  N   SER B 129      38.164  22.228  40.309  1.00 29.61           N  
ATOM   2321  CA  SER B 129      37.666  22.705  41.589  1.00 30.01           C  
ATOM   2322  C   SER B 129      36.361  21.993  42.001  1.00 30.37           C  
ATOM   2323  O   SER B 129      36.242  21.543  43.138  1.00 31.09           O  
ATOM   2324  CB  SER B 129      37.484  24.230  41.579  1.00 30.01           C  
ATOM   2325  OG  SER B 129      36.828  24.692  42.754  1.00 29.94           O  
ATOM   2326  N   VAL B 130      35.407  21.867  41.082  1.00 30.23           N  
ATOM   2327  CA  VAL B 130      34.141  21.192  41.371  1.00 30.20           C  
ATOM   2328  C   VAL B 130      34.332  19.704  41.685  1.00 30.46           C  
ATOM   2329  O   VAL B 130      33.817  19.206  42.697  1.00 30.53           O  
ATOM   2330  CB  VAL B 130      33.100  21.413  40.230  1.00 30.41           C  
ATOM   2331  CG1 VAL B 130      32.022  20.338  40.226  1.00 29.61           C  
ATOM   2332  CG2 VAL B 130      32.470  22.795  40.364  1.00 30.20           C  
ATOM   2333  N   THR B 131      35.079  18.988  40.849  1.00 30.50           N  
ATOM   2334  CA  THR B 131      35.277  17.574  41.133  1.00 30.69           C  
ATOM   2335  C   THR B 131      36.020  17.381  42.450  1.00 30.67           C  
ATOM   2336  O   THR B 131      35.805  16.372  43.127  1.00 30.67           O  
ATOM   2337  CB  THR B 131      35.971  16.761  40.005  1.00 30.67           C  
ATOM   2338  OG1 THR B 131      37.320  16.459  40.377  1.00 31.71           O  
ATOM   2339  CG2 THR B 131      35.929  17.471  38.672  1.00 30.56           C  
ATOM   2340  N   ASP B 132      36.874  18.345  42.817  1.00 30.58           N  
ATOM   2341  CA  ASP B 132      37.604  18.281  44.095  1.00 30.63           C  
ATOM   2342  C   ASP B 132      36.647  18.155  45.285  1.00 30.67           C  
ATOM   2343  O   ASP B 132      36.963  17.499  46.268  1.00 30.52           O  
ATOM   2344  CB  ASP B 132      38.509  19.505  44.284  1.00 30.65           C  
ATOM   2345  CG  ASP B 132      39.880  19.352  43.630  1.00 30.53           C  
ATOM   2346  OD1 ASP B 132      40.150  18.339  42.955  1.00 30.85           O  
ATOM   2347  OD2 ASP B 132      40.701  20.271  43.791  1.00 30.62           O  
ATOM   2348  N   PHE B 133      35.479  18.785  45.170  1.00 30.82           N  
ATOM   2349  CA  PHE B 133      34.445  18.776  46.204  1.00 31.08           C  
ATOM   2350  C   PHE B 133      33.484  17.595  46.049  1.00 30.86           C  
ATOM   2351  O   PHE B 133      32.478  17.507  46.756  1.00 30.88           O  
ATOM   2352  CB  PHE B 133      33.641  20.085  46.142  1.00 31.48           C  
ATOM   2353  CG  PHE B 133      34.341  21.260  46.759  1.00 32.27           C  
ATOM   2354  CD1 PHE B 133      34.048  21.650  48.062  1.00 33.38           C  
ATOM   2355  CD2 PHE B 133      35.299  21.974  46.045  1.00 32.82           C  
ATOM   2356  CE1 PHE B 133      34.702  22.737  48.646  1.00 34.69           C  
ATOM   2357  CE2 PHE B 133      35.950  23.068  46.612  1.00 33.36           C  
ATOM   2358  CZ  PHE B 133      35.658  23.446  47.916  1.00 33.72           C  
ATOM   2359  N   GLY B 134      33.786  16.706  45.110  1.00 30.43           N  
ATOM   2360  CA  GLY B 134      32.888  15.612  44.776  1.00 30.19           C  
ATOM   2361  C   GLY B 134      31.660  16.077  44.019  1.00 30.12           C  
ATOM   2362  O   GLY B 134      30.569  15.547  44.221  1.00 29.96           O  
ATOM   2363  N   GLY B 135      31.831  17.082  43.163  1.00 30.19           N  
ATOM   2364  CA  GLY B 135      30.756  17.517  42.277  1.00 30.44           C  
ATOM   2365  C   GLY B 135      30.103  18.833  42.640  1.00 30.74           C  
ATOM   2366  O   GLY B 135      30.381  19.418  43.690  1.00 31.14           O  
ATOM   2367  N   LEU B 136      29.208  19.283  41.767  1.00 30.67           N  
ATOM   2368  CA  LEU B 136      28.629  20.619  41.849  1.00 30.59           C  
ATOM   2369  C   LEU B 136      27.801  20.879  43.104  1.00 30.71           C  
ATOM   2370  O   LEU B 136      28.017  21.889  43.786  1.00 30.75           O  
ATOM   2371  CB  LEU B 136      27.810  20.939  40.590  1.00 30.21           C  
ATOM   2372  CG  LEU B 136      27.235  22.357  40.469  1.00 30.46           C  
ATOM   2373  CD1 LEU B 136      28.339  23.421  40.515  1.00 29.94           C  
ATOM   2374  CD2 LEU B 136      26.401  22.503  39.197  1.00 30.70           C  
ATOM   2375  N   GLU B 137      26.854  19.985  43.397  1.00 30.65           N  
ATOM   2376  CA  GLU B 137      25.934  20.186  44.518  1.00 30.70           C  
ATOM   2377  C   GLU B 137      26.690  20.318  45.837  1.00 30.12           C  
ATOM   2378  O   GLU B 137      26.291  21.092  46.690  1.00 30.16           O  
ATOM   2379  CB  GLU B 137      24.873  19.075  44.615  1.00 31.15           C  
ATOM   2380  CG  GLU B 137      24.078  18.746  43.310  1.00 33.09           C  
ATOM   2381  CD  GLU B 137      23.097  19.832  42.855  1.00 34.36           C  
ATOM   2382  OE1 GLU B 137      23.462  21.026  42.859  1.00 36.12           O  
ATOM   2383  OE2 GLU B 137      21.965  19.487  42.453  1.00 34.65           O  
ATOM   2384  N   ASN B 138      27.783  19.576  45.998  1.00 29.67           N  
ATOM   2385  CA  ASN B 138      28.629  19.726  47.176  1.00 29.43           C  
ATOM   2386  C   ASN B 138      29.264  21.107  47.233  1.00 30.00           C  
ATOM   2387  O   ASN B 138      29.231  21.788  48.277  1.00 29.74           O  
ATOM   2388  CB  ASN B 138      29.714  18.664  47.204  1.00 29.21           C  
ATOM   2389  CG  ASN B 138      29.222  17.338  47.744  1.00 28.26           C  
ATOM   2390  OD1 ASN B 138      28.032  17.160  48.024  1.00 27.38           O  
ATOM   2391  ND2 ASN B 138      30.138  16.394  47.890  1.00 25.39           N  
ATOM   2392  N   TYR B 139      29.833  21.514  46.096  1.00 30.45           N  
ATOM   2393  CA  TYR B 139      30.413  22.837  45.955  1.00 30.79           C  
ATOM   2394  C   TYR B 139      29.374  23.890  46.317  1.00 31.19           C  
ATOM   2395  O   TYR B 139      29.619  24.690  47.216  1.00 31.52           O  
ATOM   2396  CB  TYR B 139      30.955  23.060  44.547  1.00 30.79           C  
ATOM   2397  CG  TYR B 139      31.757  24.337  44.405  1.00 30.80           C  
ATOM   2398  CD1 TYR B 139      33.146  24.312  44.448  1.00 29.73           C  
ATOM   2399  CD2 TYR B 139      31.120  25.571  44.222  1.00 31.46           C  
ATOM   2400  CE1 TYR B 139      33.882  25.469  44.317  1.00 30.28           C  
ATOM   2401  CE2 TYR B 139      31.849  26.735  44.095  1.00 31.24           C  
ATOM   2402  CZ  TYR B 139      33.230  26.676  44.140  1.00 30.83           C  
ATOM   2403  OH  TYR B 139      33.960  27.831  44.006  1.00 31.42           O  
ATOM   2404  N   LYS B 140      28.221  23.866  45.642  1.00 31.27           N  
ATOM   2405  CA  LYS B 140      27.097  24.747  45.975  1.00 31.69           C  
ATOM   2406  C   LYS B 140      26.830  24.787  47.479  1.00 31.95           C  
ATOM   2407  O   LYS B 140      26.799  25.864  48.091  1.00 31.86           O  
ATOM   2408  CB  LYS B 140      25.820  24.321  45.234  1.00 31.63           C  
ATOM   2409  CG  LYS B 140      25.777  24.738  43.770  1.00 32.20           C  
ATOM   2410  CD  LYS B 140      24.483  24.311  43.094  1.00 31.95           C  
ATOM   2411  CE  LYS B 140      24.445  24.818  41.663  1.00 32.70           C  
ATOM   2412  NZ  LYS B 140      23.286  24.262  40.909  1.00 34.64           N  
ATOM   2413  N   GLU B 141      26.664  23.603  48.065  1.00 32.32           N  
ATOM   2414  CA  GLU B 141      26.338  23.455  49.480  1.00 32.80           C  
ATOM   2415  C   GLU B 141      27.366  24.159  50.383  1.00 32.59           C  
ATOM   2416  O   GLU B 141      26.996  24.881  51.294  1.00 32.28           O  
ATOM   2417  CB  GLU B 141      26.218  21.965  49.830  1.00 32.80           C  
ATOM   2418  CG  GLU B 141      25.414  21.674  51.083  1.00 34.77           C  
ATOM   2419  CD  GLU B 141      23.897  21.793  50.894  1.00 36.31           C  
ATOM   2420  OE1 GLU B 141      23.416  21.836  49.740  1.00 36.36           O  
ATOM   2421  OE2 GLU B 141      23.183  21.835  51.921  1.00 37.28           O  
ATOM   2422  N   LEU B 142      28.651  23.969  50.094  1.00 32.81           N  
ATOM   2423  CA  LEU B 142      29.717  24.479  50.949  1.00 32.89           C  
ATOM   2424  C   LEU B 142      29.985  25.954  50.719  1.00 32.77           C  
ATOM   2425  O   LEU B 142      30.425  26.664  51.628  1.00 33.37           O  
ATOM   2426  CB  LEU B 142      31.000  23.639  50.798  1.00 33.14           C  
ATOM   2427  CG  LEU B 142      30.956  22.251  51.481  1.00 33.81           C  
ATOM   2428  CD1 LEU B 142      32.120  21.363  51.081  1.00 34.95           C  
ATOM   2429  CD2 LEU B 142      30.865  22.325  53.014  1.00 32.72           C  
ATOM   2430  N   THR B 143      29.691  26.425  49.516  1.00 32.32           N  
ATOM   2431  CA  THR B 143      29.881  27.833  49.184  1.00 31.91           C  
ATOM   2432  C   THR B 143      28.551  28.608  49.178  1.00 31.86           C  
ATOM   2433  O   THR B 143      28.456  29.704  48.630  1.00 31.57           O  
ATOM   2434  CB  THR B 143      30.619  27.983  47.852  1.00 31.67           C  
ATOM   2435  OG1 THR B 143      29.868  27.330  46.824  1.00 32.03           O  
ATOM   2436  CG2 THR B 143      32.001  27.354  47.939  1.00 31.05           C  
ATOM   2437  N   GLY B 144      27.534  28.017  49.803  1.00 32.02           N  
ATOM   2438  CA  GLY B 144      26.230  28.646  49.995  1.00 31.89           C  
ATOM   2439  C   GLY B 144      25.409  28.871  48.742  1.00 31.87           C  
ATOM   2440  O   GLY B 144      24.809  29.929  48.589  1.00 31.94           O  
ATOM   2441  N   GLY B 145      25.387  27.884  47.845  1.00 31.80           N  
ATOM   2442  CA  GLY B 145      24.609  27.958  46.599  1.00 31.57           C  
ATOM   2443  C   GLY B 145      25.357  28.391  45.344  1.00 31.64           C  
ATOM   2444  O   GLY B 145      24.779  28.422  44.263  1.00 31.73           O  
ATOM   2445  N   ALA B 146      26.641  28.715  45.476  1.00 31.78           N  
ATOM   2446  CA  ALA B 146      27.426  29.314  44.379  1.00 32.01           C  
ATOM   2447  C   ALA B 146      27.735  28.374  43.212  1.00 31.79           C  
ATOM   2448  O   ALA B 146      28.488  27.416  43.369  1.00 32.07           O  
ATOM   2449  CB  ALA B 146      28.740  29.917  44.930  1.00 32.13           C  
ATOM   2450  N   ASP B 147      27.166  28.664  42.043  1.00 31.80           N  
ATOM   2451  CA  ASP B 147      27.500  27.947  40.803  1.00 31.39           C  
ATOM   2452  C   ASP B 147      28.709  28.639  40.189  1.00 31.14           C  
ATOM   2453  O   ASP B 147      28.587  29.748  39.659  1.00 31.10           O  
ATOM   2454  CB  ASP B 147      26.315  27.950  39.830  1.00 31.41           C  
ATOM   2455  CG  ASP B 147      26.605  27.204  38.514  1.00 32.72           C  
ATOM   2456  OD1 ASP B 147      27.705  26.635  38.323  1.00 33.05           O  
ATOM   2457  OD2 ASP B 147      25.706  27.189  37.646  1.00 34.58           O  
ATOM   2458  N   PRO B 148      29.884  27.990  40.260  1.00 30.87           N  
ATOM   2459  CA  PRO B 148      31.120  28.646  39.839  1.00 30.73           C  
ATOM   2460  C   PRO B 148      31.246  28.749  38.320  1.00 30.71           C  
ATOM   2461  O   PRO B 148      32.163  29.407  37.832  1.00 30.66           O  
ATOM   2462  CB  PRO B 148      32.208  27.738  40.414  1.00 30.71           C  
ATOM   2463  CG  PRO B 148      31.580  26.392  40.453  1.00 30.95           C  
ATOM   2464  CD  PRO B 148      30.116  26.605  40.710  1.00 30.76           C  
ATOM   2465  N   PHE B 149      30.329  28.099  37.596  1.00 30.75           N  
ATOM   2466  CA  PHE B 149      30.282  28.130  36.127  1.00 30.79           C  
ATOM   2467  C   PHE B 149      29.442  29.299  35.623  1.00 31.21           C  
ATOM   2468  O   PHE B 149      29.417  29.596  34.418  1.00 31.14           O  
ATOM   2469  CB  PHE B 149      29.701  26.822  35.561  1.00 30.17           C  
ATOM   2470  CG  PHE B 149      30.625  25.643  35.661  1.00 29.43           C  
ATOM   2471  CD1 PHE B 149      31.726  25.526  34.818  1.00 28.73           C  
ATOM   2472  CD2 PHE B 149      30.385  24.633  36.585  1.00 28.64           C  
ATOM   2473  CE1 PHE B 149      32.577  24.431  34.907  1.00 28.23           C  
ATOM   2474  CE2 PHE B 149      31.232  23.535  36.681  1.00 27.39           C  
ATOM   2475  CZ  PHE B 149      32.326  23.433  35.843  1.00 28.35           C  
ATOM   2476  N   ALA B 150      28.744  29.947  36.549  1.00 31.81           N  
ATOM   2477  CA  ALA B 150      27.844  31.047  36.211  1.00 32.43           C  
ATOM   2478  C   ALA B 150      28.597  32.294  35.742  1.00 32.80           C  
ATOM   2479  O   ALA B 150      29.763  32.497  36.094  1.00 32.41           O  
ATOM   2480  CB  ALA B 150      26.945  31.379  37.404  1.00 32.12           C  
ATOM   2481  N   LEU B 151      27.921  33.114  34.940  1.00 33.62           N  
ATOM   2482  CA  LEU B 151      28.431  34.425  34.570  1.00 34.71           C  
ATOM   2483  C   LEU B 151      28.400  35.322  35.803  1.00 35.59           C  
ATOM   2484  O   LEU B 151      27.394  35.394  36.512  1.00 35.34           O  
ATOM   2485  CB  LEU B 151      27.604  35.045  33.434  1.00 34.71           C  
ATOM   2486  CG  LEU B 151      27.170  34.192  32.228  1.00 34.53           C  
ATOM   2487  CD1 LEU B 151      26.125  34.943  31.380  1.00 34.42           C  
ATOM   2488  CD2 LEU B 151      28.352  33.725  31.377  1.00 33.17           C  
HETATM 2489  N   MSE B 152      29.524  35.982  36.060  1.00 36.98           N  
HETATM 2490  CA  MSE B 152      29.677  36.849  37.223  1.00 38.36           C  
HETATM 2491  C   MSE B 152      29.327  38.297  36.915  1.00 37.79           C  
HETATM 2492  O   MSE B 152      29.603  38.796  35.814  1.00 37.92           O  
HETATM 2493  CB  MSE B 152      31.103  36.766  37.778  1.00 39.17           C  
HETATM 2494  CG  MSE B 152      31.542  35.365  38.214  1.00 44.43           C  
HETATM 2495 SE   MSE B 152      30.186  34.348  39.245  1.00 58.74          SE  
HETATM 2496  CE  MSE B 152      29.769  35.652  40.669  1.00 54.99           C  
ATOM   2497  N   THR B 153      28.681  38.947  37.880  1.00 37.25           N  
ATOM   2498  CA  THR B 153      28.610  40.401  37.925  1.00 36.87           C  
ATOM   2499  C   THR B 153      29.616  40.791  38.998  1.00 36.22           C  
ATOM   2500  O   THR B 153      29.411  40.503  40.169  1.00 36.12           O  
ATOM   2501  CB  THR B 153      27.172  40.945  38.238  1.00 37.09           C  
ATOM   2502  OG1 THR B 153      26.763  40.556  39.552  1.00 36.41           O  
ATOM   2503  CG2 THR B 153      26.151  40.416  37.232  1.00 37.96           C  
ATOM   2504  N   PRO B 154      30.727  41.430  38.601  1.00 36.11           N  
ATOM   2505  CA  PRO B 154      31.863  41.565  39.530  1.00 35.88           C  
ATOM   2506  C   PRO B 154      31.629  42.570  40.665  1.00 35.54           C  
ATOM   2507  O   PRO B 154      31.095  43.658  40.423  1.00 35.75           O  
ATOM   2508  CB  PRO B 154      32.991  42.062  38.623  1.00 35.82           C  
ATOM   2509  CG  PRO B 154      32.263  42.846  37.549  1.00 35.84           C  
ATOM   2510  CD  PRO B 154      30.994  42.077  37.300  1.00 35.91           C  
ATOM   2511  N   VAL B 155      32.025  42.194  41.881  1.00 34.83           N  
ATOM   2512  CA  VAL B 155      32.149  43.126  43.007  1.00 34.06           C  
ATOM   2513  C   VAL B 155      33.088  44.272  42.620  1.00 33.95           C  
ATOM   2514  O   VAL B 155      33.849  44.156  41.649  1.00 33.68           O  
ATOM   2515  CB  VAL B 155      32.677  42.421  44.299  1.00 34.28           C  
ATOM   2516  CG1 VAL B 155      31.718  41.324  44.753  1.00 32.87           C  
ATOM   2517  CG2 VAL B 155      34.093  41.845  44.090  1.00 34.22           C  
HETATM 2518  N   OCS B 156      33.037  45.369  43.378  1.00 33.72           N  
HETATM 2519  CA  OCS B 156      33.779  46.584  43.034  1.00 33.44           C  
HETATM 2520  CB  OCS B 156      33.456  47.727  43.997  1.00 33.65           C  
HETATM 2521  SG  OCS B 156      33.767  47.393  45.739  1.00 32.53           S  
HETATM 2522  C   OCS B 156      35.291  46.402  42.914  1.00 33.60           C  
HETATM 2523  O   OCS B 156      35.900  46.997  42.029  1.00 34.00           O  
HETATM 2524  OD1 OCS B 156      35.145  47.608  46.048  1.00 31.18           O  
HETATM 2525  OD2 OCS B 156      32.950  48.213  46.556  1.00 32.62           O  
HETATM 2526  OD3 OCS B 156      33.403  46.054  45.946  1.00 34.90           O  
ATOM   2527  N   GLY B 157      35.887  45.591  43.788  1.00 33.21           N  
ATOM   2528  CA  GLY B 157      37.304  45.269  43.691  1.00 33.28           C  
ATOM   2529  C   GLY B 157      37.617  44.585  42.373  1.00 33.72           C  
ATOM   2530  O   GLY B 157      38.523  44.992  41.651  1.00 33.78           O  
ATOM   2531  N   LEU B 158      36.836  43.558  42.048  1.00 34.11           N  
ATOM   2532  CA  LEU B 158      37.033  42.771  40.836  1.00 34.48           C  
ATOM   2533  C   LEU B 158      36.901  43.606  39.541  1.00 35.12           C  
ATOM   2534  O   LEU B 158      37.796  43.568  38.689  1.00 35.67           O  
ATOM   2535  CB  LEU B 158      36.096  41.551  40.844  1.00 34.42           C  
ATOM   2536  CG  LEU B 158      36.351  40.338  39.930  1.00 34.65           C  
ATOM   2537  CD1 LEU B 158      37.729  39.704  40.152  1.00 34.43           C  
ATOM   2538  CD2 LEU B 158      35.246  39.283  40.086  1.00 33.99           C  
ATOM   2539  N   SER B 159      35.811  44.360  39.393  1.00 35.23           N  
ATOM   2540  CA  SER B 159      35.634  45.222  38.216  1.00 35.53           C  
ATOM   2541  C   SER B 159      36.759  46.249  38.094  1.00 35.88           C  
ATOM   2542  O   SER B 159      37.363  46.392  37.033  1.00 36.37           O  
ATOM   2543  CB  SER B 159      34.271  45.936  38.230  1.00 35.32           C  
ATOM   2544  OG  SER B 159      34.106  46.727  39.402  1.00 35.07           O  
ATOM   2545  N   ALA B 160      37.039  46.959  39.182  1.00 36.23           N  
ATOM   2546  CA  ALA B 160      38.087  47.980  39.182  1.00 36.32           C  
ATOM   2547  C   ALA B 160      39.455  47.409  38.766  1.00 36.29           C  
ATOM   2548  O   ALA B 160      40.073  47.912  37.837  1.00 36.56           O  
ATOM   2549  CB  ALA B 160      38.166  48.673  40.540  1.00 36.01           C  
ATOM   2550  N   ASN B 161      39.905  46.348  39.429  1.00 36.10           N  
ATOM   2551  CA  ASN B 161      41.196  45.748  39.112  1.00 36.18           C  
ATOM   2552  C   ASN B 161      41.308  45.254  37.673  1.00 36.34           C  
ATOM   2553  O   ASN B 161      42.367  45.370  37.045  1.00 36.68           O  
ATOM   2554  CB  ASN B 161      41.504  44.593  40.059  1.00 36.15           C  
ATOM   2555  CG  ASN B 161      41.782  45.053  41.484  1.00 36.24           C  
ATOM   2556  OD1 ASN B 161      42.350  46.120  41.707  1.00 35.99           O  
ATOM   2557  ND2 ASN B 161      41.396  44.229  42.455  1.00 34.98           N  
ATOM   2558  N   ASN B 162      40.222  44.690  37.154  1.00 36.19           N  
ATOM   2559  CA  ASN B 162      40.244  44.113  35.818  1.00 35.61           C  
ATOM   2560  C   ASN B 162      40.144  45.171  34.724  1.00 35.40           C  
ATOM   2561  O   ASN B 162      40.866  45.098  33.729  1.00 35.46           O  
ATOM   2562  CB  ASN B 162      39.156  43.058  35.670  1.00 35.71           C  
ATOM   2563  CG  ASN B 162      39.559  41.722  36.263  1.00 36.51           C  
ATOM   2564  OD1 ASN B 162      40.577  41.146  35.885  1.00 36.94           O  
ATOM   2565  ND2 ASN B 162      38.750  41.212  37.186  1.00 37.94           N  
ATOM   2566  N   ILE B 163      39.264  46.158  34.910  1.00 34.90           N  
ATOM   2567  CA  ILE B 163      39.160  47.252  33.949  1.00 34.52           C  
ATOM   2568  C   ILE B 163      40.448  48.100  33.944  1.00 34.72           C  
ATOM   2569  O   ILE B 163      40.888  48.551  32.887  1.00 34.82           O  
ATOM   2570  CB  ILE B 163      37.841  48.066  34.103  1.00 34.47           C  
ATOM   2571  CG1 ILE B 163      36.654  47.189  33.682  1.00 34.57           C  
ATOM   2572  CG2 ILE B 163      37.872  49.306  33.251  1.00 33.05           C  
ATOM   2573  CD1 ILE B 163      35.275  47.788  33.884  1.00 33.91           C  
ATOM   2574  N   PHE B 164      41.070  48.277  35.108  1.00 34.77           N  
ATOM   2575  CA  PHE B 164      42.383  48.930  35.190  1.00 34.73           C  
ATOM   2576  C   PHE B 164      43.409  48.227  34.291  1.00 35.00           C  
ATOM   2577  O   PHE B 164      44.066  48.878  33.467  1.00 34.67           O  
ATOM   2578  CB  PHE B 164      42.893  48.982  36.633  1.00 34.09           C  
ATOM   2579  CG  PHE B 164      44.145  49.808  36.809  1.00 33.92           C  
ATOM   2580  CD1 PHE B 164      44.064  51.160  37.141  1.00 32.83           C  
ATOM   2581  CD2 PHE B 164      45.406  49.234  36.656  1.00 33.15           C  
ATOM   2582  CE1 PHE B 164      45.212  51.925  37.308  1.00 31.92           C  
ATOM   2583  CE2 PHE B 164      46.560  49.992  36.821  1.00 32.07           C  
ATOM   2584  CZ  PHE B 164      46.461  51.338  37.149  1.00 32.61           C  
ATOM   2585  N   LYS B 165      43.521  46.904  34.454  1.00 35.36           N  
ATOM   2586  CA  LYS B 165      44.449  46.064  33.678  1.00 35.67           C  
ATOM   2587  C   LYS B 165      44.168  46.147  32.190  1.00 36.08           C  
ATOM   2588  O   LYS B 165      45.087  46.307  31.390  1.00 35.89           O  
ATOM   2589  CB  LYS B 165      44.343  44.598  34.107  1.00 35.61           C  
ATOM   2590  CG  LYS B 165      45.012  44.259  35.422  1.00 35.15           C  
ATOM   2591  CD  LYS B 165      44.785  42.797  35.804  1.00 35.32           C  
ATOM   2592  CE  LYS B 165      45.449  41.813  34.833  1.00 33.92           C  
ATOM   2593  NZ  LYS B 165      46.926  41.709  35.008  1.00 32.67           N  
ATOM   2594  N   LEU B 166      42.886  46.039  31.839  1.00 36.82           N  
ATOM   2595  CA  LEU B 166      42.421  46.099  30.455  1.00 37.50           C  
ATOM   2596  C   LEU B 166      42.807  47.397  29.756  1.00 38.37           C  
ATOM   2597  O   LEU B 166      43.109  47.393  28.570  1.00 38.62           O  
ATOM   2598  CB  LEU B 166      40.900  45.914  30.397  1.00 37.37           C  
ATOM   2599  CG  LEU B 166      40.215  45.707  29.040  1.00 37.03           C  
ATOM   2600  CD1 LEU B 166      40.608  44.396  28.394  1.00 36.27           C  
ATOM   2601  CD2 LEU B 166      38.723  45.731  29.234  1.00 37.20           C  
HETATM 2602  N   MSE B 167      42.804  48.500  30.496  1.00 39.25           N  
HETATM 2603  CA  MSE B 167      43.066  49.809  29.917  1.00 40.95           C  
HETATM 2604  C   MSE B 167      44.549  50.084  29.821  1.00 39.30           C  
HETATM 2605  O   MSE B 167      45.019  50.619  28.826  1.00 39.67           O  
HETATM 2606  CB  MSE B 167      42.356  50.912  30.713  1.00 40.43           C  
HETATM 2607  CG  MSE B 167      40.838  50.882  30.572  1.00 42.21           C  
HETATM 2608 SE   MSE B 167      39.868  52.288  31.541  1.00 48.11          SE  
HETATM 2609  CE  MSE B 167      40.439  53.837  30.535  1.00 45.79           C  
ATOM   2610  N   THR B 168      45.288  49.696  30.846  1.00 38.64           N  
ATOM   2611  CA  THR B 168      46.689  50.086  30.970  1.00 37.80           C  
ATOM   2612  C   THR B 168      47.658  49.152  30.254  1.00 37.53           C  
ATOM   2613  O   THR B 168      48.676  49.605  29.727  1.00 37.45           O  
ATOM   2614  CB  THR B 168      47.116  50.213  32.451  1.00 37.65           C  
ATOM   2615  OG1 THR B 168      46.806  49.003  33.152  1.00 37.17           O  
ATOM   2616  CG2 THR B 168      46.397  51.356  33.114  1.00 37.17           C  
ATOM   2617  N   GLU B 169      47.330  47.860  30.227  1.00 37.26           N  
ATOM   2618  CA  GLU B 169      48.274  46.810  29.798  1.00 36.96           C  
ATOM   2619  C   GLU B 169      48.650  46.811  28.322  1.00 36.47           C  
ATOM   2620  O   GLU B 169      47.803  46.950  27.446  1.00 36.26           O  
ATOM   2621  CB  GLU B 169      47.793  45.420  30.231  1.00 36.88           C  
ATOM   2622  CG  GLU B 169      48.113  45.128  31.691  1.00 37.13           C  
ATOM   2623  CD  GLU B 169      47.604  43.787  32.192  1.00 37.13           C  
ATOM   2624  OE1 GLU B 169      47.129  42.954  31.390  1.00 37.65           O  
ATOM   2625  OE2 GLU B 169      47.700  43.568  33.413  1.00 38.07           O  
ATOM   2626  N   LYS B 170      49.946  46.667  28.077  1.00 36.34           N  
ATOM   2627  CA  LYS B 170      50.488  46.548  26.732  1.00 36.25           C  
ATOM   2628  C   LYS B 170      51.002  45.117  26.545  1.00 36.07           C  
ATOM   2629  O   LYS B 170      50.277  44.236  26.071  1.00 35.98           O  
ATOM   2630  CB  LYS B 170      51.618  47.558  26.524  1.00 36.09           C  
ATOM   2631  CG  LYS B 170      52.009  47.764  25.071  1.00 36.71           C  
ATOM   2632  CD  LYS B 170      51.461  49.079  24.521  1.00 37.05           C  
ATOM   2633  CE  LYS B 170      51.679  49.190  23.020  1.00 37.18           C  
ATOM   2634  NZ  LYS B 170      52.979  48.616  22.588  1.00 37.77           N  
ATOM   2635  N   ASP B 171      52.247  44.892  26.953  1.00 35.78           N  
ATOM   2636  CA  ASP B 171      52.892  43.597  26.810  1.00 35.43           C  
ATOM   2637  C   ASP B 171      53.451  43.099  28.148  1.00 35.46           C  
ATOM   2638  O   ASP B 171      54.090  42.046  28.216  1.00 35.63           O  
ATOM   2639  CB  ASP B 171      53.983  43.678  25.738  1.00 35.16           C  
ATOM   2640  CG  ASP B 171      55.017  44.761  26.023  1.00 34.80           C  
ATOM   2641  OD1 ASP B 171      54.749  45.683  26.826  1.00 33.12           O  
ATOM   2642  OD2 ASP B 171      56.113  44.686  25.429  1.00 34.69           O  
ATOM   2643  N   VAL B 172      53.190  43.860  29.208  1.00 35.48           N  
ATOM   2644  CA  VAL B 172      53.694  43.557  30.546  1.00 35.44           C  
ATOM   2645  C   VAL B 172      52.540  43.541  31.548  1.00 35.53           C  
ATOM   2646  O   VAL B 172      51.745  44.483  31.599  1.00 35.35           O  
ATOM   2647  CB  VAL B 172      54.776  44.582  30.999  1.00 35.48           C  
ATOM   2648  CG1 VAL B 172      55.163  44.362  32.460  1.00 35.68           C  
ATOM   2649  CG2 VAL B 172      56.008  44.517  30.098  1.00 34.90           C  
ATOM   2650  N   PRO B 173      52.437  42.462  32.344  1.00 35.73           N  
ATOM   2651  CA  PRO B 173      51.421  42.405  33.392  1.00 35.84           C  
ATOM   2652  C   PRO B 173      51.461  43.651  34.278  1.00 36.04           C  
ATOM   2653  O   PRO B 173      52.534  44.094  34.686  1.00 36.20           O  
ATOM   2654  CB  PRO B 173      51.830  41.176  34.205  1.00 35.71           C  
ATOM   2655  CG  PRO B 173      52.545  40.316  33.241  1.00 35.82           C  
ATOM   2656  CD  PRO B 173      53.263  41.242  32.308  1.00 35.60           C  
ATOM   2657  N   ILE B 174      50.296  44.227  34.539  1.00 36.34           N  
ATOM   2658  CA  ILE B 174      50.185  45.322  35.493  1.00 36.60           C  
ATOM   2659  C   ILE B 174      49.249  44.881  36.624  1.00 37.03           C  
ATOM   2660  O   ILE B 174      48.031  44.817  36.443  1.00 36.89           O  
ATOM   2661  CB  ILE B 174      49.709  46.645  34.817  1.00 36.57           C  
ATOM   2662  CG1 ILE B 174      50.779  47.166  33.850  1.00 36.10           C  
ATOM   2663  CG2 ILE B 174      49.375  47.714  35.867  1.00 36.74           C  
ATOM   2664  CD1 ILE B 174      50.328  48.296  32.951  1.00 36.13           C  
ATOM   2665  N   ASP B 175      49.838  44.553  37.773  1.00 37.31           N  
ATOM   2666  CA  ASP B 175      49.084  44.123  38.950  1.00 37.96           C  
ATOM   2667  C   ASP B 175      48.596  45.337  39.762  1.00 38.13           C  
ATOM   2668  O   ASP B 175      49.386  45.981  40.456  1.00 38.31           O  
ATOM   2669  CB  ASP B 175      49.940  43.187  39.821  1.00 37.91           C  
ATOM   2670  CG  ASP B 175      49.131  42.478  40.903  1.00 38.54           C  
ATOM   2671  OD1 ASP B 175      47.944  42.820  41.100  1.00 39.34           O  
ATOM   2672  OD2 ASP B 175      49.681  41.566  41.560  1.00 39.03           O  
ATOM   2673  N   PRO B 176      47.289  45.640  39.689  1.00 38.38           N  
ATOM   2674  CA  PRO B 176      46.779  46.865  40.313  1.00 38.60           C  
ATOM   2675  C   PRO B 176      46.759  46.796  41.828  1.00 38.83           C  
ATOM   2676  O   PRO B 176      46.881  47.819  42.492  1.00 38.96           O  
ATOM   2677  CB  PRO B 176      45.362  47.014  39.747  1.00 38.68           C  
ATOM   2678  CG  PRO B 176      45.025  45.720  39.082  1.00 38.79           C  
ATOM   2679  CD  PRO B 176      46.236  44.847  39.024  1.00 38.70           C  
ATOM   2680  N   THR B 177      46.624  45.591  42.364  1.00 39.33           N  
ATOM   2681  CA  THR B 177      46.706  45.366  43.805  1.00 39.70           C  
ATOM   2682  C   THR B 177      48.160  45.346  44.293  1.00 40.04           C  
ATOM   2683  O   THR B 177      48.451  44.791  45.348  1.00 40.33           O  
ATOM   2684  CB  THR B 177      46.050  44.024  44.202  1.00 39.77           C  
ATOM   2685  OG1 THR B 177      46.879  42.940  43.763  1.00 39.47           O  
ATOM   2686  CG2 THR B 177      44.651  43.891  43.596  1.00 39.00           C  
ATOM   2687  N   SER B 178      49.065  45.956  43.532  1.00 40.42           N  
ATOM   2688  CA  SER B 178      50.500  45.854  43.801  1.00 40.60           C  
ATOM   2689  C   SER B 178      51.280  47.119  43.424  1.00 40.67           C  
ATOM   2690  O   SER B 178      52.483  47.210  43.671  1.00 40.53           O  
ATOM   2691  CB  SER B 178      51.063  44.613  43.089  1.00 40.69           C  
ATOM   2692  OG  SER B 178      52.325  44.851  42.488  1.00 41.19           O  
ATOM   2693  N   ILE B 179      50.591  48.091  42.833  1.00 40.88           N  
ATOM   2694  CA  ILE B 179      51.221  49.354  42.437  1.00 41.20           C  
ATOM   2695  C   ILE B 179      51.242  50.336  43.608  1.00 41.51           C  
ATOM   2696  O   ILE B 179      50.446  50.208  44.536  1.00 41.53           O  
ATOM   2697  CB  ILE B 179      50.539  49.986  41.181  1.00 41.18           C  
ATOM   2698  CG1 ILE B 179      49.045  50.244  41.429  1.00 41.14           C  
ATOM   2699  CG2 ILE B 179      50.749  49.095  39.954  1.00 40.87           C  
ATOM   2700  CD1 ILE B 179      48.360  51.060  40.352  1.00 40.98           C  
ATOM   2701  N   GLU B 180      52.160  51.300  43.567  1.00 41.92           N  
ATOM   2702  CA  GLU B 180      52.279  52.311  44.623  1.00 42.46           C  
ATOM   2703  C   GLU B 180      51.083  53.264  44.599  1.00 42.56           C  
ATOM   2704  O   GLU B 180      50.839  53.943  43.593  1.00 42.71           O  
ATOM   2705  CB  GLU B 180      53.584  53.101  44.477  1.00 42.26           C  
ATOM   2706  CG  GLU B 180      54.856  52.270  44.617  1.00 42.65           C  
ATOM   2707  CD  GLU B 180      56.130  53.063  44.321  1.00 43.17           C  
ATOM   2708  OE1 GLU B 180      57.232  52.522  44.565  1.00 43.99           O  
ATOM   2709  OE2 GLU B 180      56.044  54.222  43.847  1.00 43.93           O  
ATOM   2710  N   TYR B 181      50.334  53.294  45.702  1.00 42.57           N  
ATOM   2711  CA  TYR B 181      49.178  54.183  45.840  1.00 42.55           C  
ATOM   2712  C   TYR B 181      49.564  55.419  46.637  1.00 42.68           C  
ATOM   2713  O   TYR B 181      50.215  55.313  47.681  1.00 42.67           O  
ATOM   2714  CB  TYR B 181      48.032  53.475  46.567  1.00 42.53           C  
ATOM   2715  CG  TYR B 181      47.238  52.478  45.752  1.00 42.40           C  
ATOM   2716  CD1 TYR B 181      45.930  52.761  45.359  1.00 42.80           C  
ATOM   2717  CD2 TYR B 181      47.780  51.241  45.397  1.00 42.44           C  
ATOM   2718  CE1 TYR B 181      45.184  51.841  44.613  1.00 42.76           C  
ATOM   2719  CE2 TYR B 181      47.048  50.316  44.651  1.00 42.49           C  
ATOM   2720  CZ  TYR B 181      45.750  50.623  44.263  1.00 42.62           C  
ATOM   2721  OH  TYR B 181      45.022  49.721  43.525  1.00 42.33           O  
ATOM   2722  N   LEU B 182      49.159  56.590  46.153  1.00 42.84           N  
ATOM   2723  CA  LEU B 182      49.386  57.840  46.884  1.00 42.93           C  
ATOM   2724  C   LEU B 182      48.378  57.991  48.018  1.00 42.93           C  
ATOM   2725  O   LEU B 182      47.256  57.481  47.935  1.00 42.94           O  
ATOM   2726  CB  LEU B 182      49.287  59.046  45.945  1.00 43.01           C  
ATOM   2727  CG  LEU B 182      50.298  59.176  44.802  1.00 43.18           C  
ATOM   2728  CD1 LEU B 182      49.720  60.042  43.704  1.00 43.05           C  
ATOM   2729  CD2 LEU B 182      51.631  59.741  45.295  1.00 43.72           C  
ATOM   2730  N   GLU B 183      48.789  58.685  49.075  1.00 42.91           N  
ATOM   2731  CA  GLU B 183      47.901  59.023  50.185  1.00 42.94           C  
ATOM   2732  C   GLU B 183      48.275  60.387  50.744  1.00 42.81           C  
ATOM   2733  O   GLU B 183      47.430  61.098  51.289  1.00 42.83           O  
ATOM   2734  CB  GLU B 183      47.968  57.962  51.288  1.00 42.95           C  
ATOM   2735  CG  GLU B 183      46.797  58.023  52.272  1.00 43.18           C  
ATOM   2736  CD  GLU B 183      46.959  57.097  53.473  1.00 43.19           C  
ATOM   2737  OE1 GLU B 183      47.449  55.956  53.302  1.00 43.82           O  
ATOM   2738  OE2 GLU B 183      46.582  57.512  54.592  1.00 42.74           O  
ATOM   2739  N   ASN B 184      49.547  60.746  50.579  1.00 42.76           N  
ATOM   2740  CA  ASN B 184      50.113  61.986  51.116  1.00 42.71           C  
ATOM   2741  C   ASN B 184      49.792  63.220  50.271  1.00 42.37           C  
ATOM   2742  O   ASN B 184      50.322  64.310  50.522  1.00 42.20           O  
ATOM   2743  CB  ASN B 184      51.639  61.855  51.259  1.00 42.95           C  
ATOM   2744  CG  ASN B 184      52.059  60.746  52.220  1.00 43.30           C  
ATOM   2745  OD1 ASN B 184      52.784  60.998  53.184  1.00 43.34           O  
ATOM   2746  ND2 ASN B 184      51.621  59.514  51.950  1.00 43.49           N  
ATOM   2747  N   THR B 185      48.929  63.042  49.273  1.00 42.02           N  
ATOM   2748  CA  THR B 185      48.619  64.106  48.320  1.00 41.61           C  
ATOM   2749  C   THR B 185      47.130  64.163  47.971  1.00 41.15           C  
ATOM   2750  O   THR B 185      46.468  63.128  47.882  1.00 41.16           O  
ATOM   2751  CB  THR B 185      49.508  63.999  47.038  1.00 41.55           C  
ATOM   2752  OG1 THR B 185      49.499  65.257  46.341  1.00 42.02           O  
ATOM   2753  CG2 THR B 185      49.012  62.911  46.102  1.00 41.48           C  
ATOM   2754  N   SER B 186      46.612  65.377  47.794  1.00 40.66           N  
ATOM   2755  CA  SER B 186      45.214  65.575  47.417  1.00 40.23           C  
ATOM   2756  C   SER B 186      45.011  65.240  45.940  1.00 40.14           C  
ATOM   2757  O   SER B 186      45.981  65.058  45.190  1.00 40.03           O  
ATOM   2758  CB  SER B 186      44.771  67.012  47.688  1.00 40.15           C  
ATOM   2759  OG  SER B 186      45.002  67.828  46.555  1.00 39.72           O  
ATOM   2760  N   PHE B 187      43.747  65.170  45.528  1.00 39.77           N  
ATOM   2761  CA  PHE B 187      43.419  64.805  44.159  1.00 39.38           C  
ATOM   2762  C   PHE B 187      43.773  65.891  43.149  1.00 39.32           C  
ATOM   2763  O   PHE B 187      44.358  65.596  42.103  1.00 39.26           O  
ATOM   2764  CB  PHE B 187      41.949  64.405  44.023  1.00 39.24           C  
ATOM   2765  CG  PHE B 187      41.586  63.920  42.650  1.00 39.06           C  
ATOM   2766  CD1 PHE B 187      42.294  62.870  42.056  1.00 38.58           C  
ATOM   2767  CD2 PHE B 187      40.547  64.507  41.945  1.00 38.80           C  
ATOM   2768  CE1 PHE B 187      41.971  62.418  40.783  1.00 38.28           C  
ATOM   2769  CE2 PHE B 187      40.211  64.053  40.667  1.00 38.92           C  
ATOM   2770  CZ  PHE B 187      40.925  63.007  40.089  1.00 38.71           C  
ATOM   2771  N   ALA B 188      43.418  67.136  43.467  1.00 39.16           N  
ATOM   2772  CA  ALA B 188      43.711  68.279  42.600  1.00 39.04           C  
ATOM   2773  C   ALA B 188      45.197  68.372  42.253  1.00 38.75           C  
ATOM   2774  O   ALA B 188      45.560  68.499  41.081  1.00 38.84           O  
ATOM   2775  CB  ALA B 188      43.238  69.575  43.246  1.00 39.21           C  
ATOM   2776  N   GLU B 189      46.044  68.291  43.276  1.00 38.28           N  
ATOM   2777  CA  GLU B 189      47.488  68.423  43.111  1.00 37.99           C  
ATOM   2778  C   GLU B 189      48.050  67.397  42.138  1.00 37.51           C  
ATOM   2779  O   GLU B 189      48.934  67.713  41.348  1.00 37.57           O  
ATOM   2780  CB  GLU B 189      48.206  68.339  44.465  1.00 37.98           C  
ATOM   2781  CG  GLU B 189      48.026  69.574  45.344  1.00 38.32           C  
ATOM   2782  CD  GLU B 189      48.461  69.347  46.786  1.00 38.32           C  
ATOM   2783  OE1 GLU B 189      49.544  69.841  47.156  1.00 38.69           O  
ATOM   2784  OE2 GLU B 189      47.725  68.676  47.545  1.00 38.28           O  
ATOM   2785  N   HIS B 190      47.526  66.178  42.184  1.00 37.10           N  
ATOM   2786  CA  HIS B 190      47.992  65.125  41.286  1.00 36.81           C  
ATOM   2787  C   HIS B 190      47.599  65.362  39.818  1.00 36.41           C  
ATOM   2788  O   HIS B 190      48.457  65.312  38.933  1.00 36.18           O  
ATOM   2789  CB  HIS B 190      47.542  63.741  41.766  1.00 36.81           C  
ATOM   2790  CG  HIS B 190      48.250  62.617  41.079  1.00 37.40           C  
ATOM   2791  ND1 HIS B 190      49.620  62.460  41.130  1.00 37.80           N  
ATOM   2792  CD2 HIS B 190      47.783  61.606  40.309  1.00 37.84           C  
ATOM   2793  CE1 HIS B 190      49.965  61.398  40.423  1.00 38.38           C  
ATOM   2794  NE2 HIS B 190      48.870  60.860  39.917  1.00 38.52           N  
ATOM   2795  N   VAL B 191      46.314  65.629  39.577  1.00 35.97           N  
ATOM   2796  CA  VAL B 191      45.799  65.920  38.236  1.00 35.65           C  
ATOM   2797  C   VAL B 191      46.618  67.023  37.555  1.00 35.30           C  
ATOM   2798  O   VAL B 191      46.871  66.963  36.347  1.00 35.19           O  
ATOM   2799  CB  VAL B 191      44.303  66.346  38.266  1.00 35.80           C  
ATOM   2800  CG1 VAL B 191      43.705  66.302  36.866  1.00 36.09           C  
ATOM   2801  CG2 VAL B 191      43.501  65.455  39.183  1.00 35.72           C  
ATOM   2802  N   ASN B 192      47.033  68.012  38.346  1.00 34.86           N  
ATOM   2803  CA  ASN B 192      47.873  69.119  37.877  1.00 34.61           C  
ATOM   2804  C   ASN B 192      49.257  68.709  37.374  1.00 34.27           C  
ATOM   2805  O   ASN B 192      49.857  69.430  36.580  1.00 34.33           O  
ATOM   2806  CB  ASN B 192      48.028  70.178  38.974  1.00 34.63           C  
ATOM   2807  CG  ASN B 192      46.905  71.204  38.971  1.00 34.70           C  
ATOM   2808  OD1 ASN B 192      45.771  70.917  38.575  1.00 33.74           O  
ATOM   2809  ND2 ASN B 192      47.221  72.415  39.422  1.00 34.80           N  
ATOM   2810  N   THR B 193      49.761  67.567  37.845  1.00 33.92           N  
ATOM   2811  CA  THR B 193      51.080  67.065  37.434  1.00 33.38           C  
ATOM   2812  C   THR B 193      51.005  66.292  36.123  1.00 32.91           C  
ATOM   2813  O   THR B 193      52.026  66.059  35.483  1.00 33.04           O  
ATOM   2814  CB  THR B 193      51.738  66.141  38.502  1.00 33.41           C  
ATOM   2815  OG1 THR B 193      51.051  64.886  38.551  1.00 33.20           O  
ATOM   2816  CG2 THR B 193      51.731  66.785  39.889  1.00 33.31           C  
ATOM   2817  N   LEU B 194      49.794  65.899  35.735  1.00 32.25           N  
ATOM   2818  CA  LEU B 194      49.587  65.047  34.563  1.00 31.64           C  
ATOM   2819  C   LEU B 194      49.564  65.831  33.255  1.00 31.27           C  
ATOM   2820  O   LEU B 194      49.044  66.949  33.198  1.00 31.18           O  
ATOM   2821  CB  LEU B 194      48.289  64.249  34.709  1.00 31.65           C  
ATOM   2822  CG  LEU B 194      48.078  63.420  35.980  1.00 31.35           C  
ATOM   2823  CD1 LEU B 194      46.811  62.613  35.845  1.00 31.29           C  
ATOM   2824  CD2 LEU B 194      49.257  62.503  36.276  1.00 31.19           C  
ATOM   2825  N   ASP B 195      50.121  65.230  32.205  1.00 30.72           N  
ATOM   2826  CA  ASP B 195      50.174  65.855  30.884  1.00 30.16           C  
ATOM   2827  C   ASP B 195      48.800  65.844  30.200  1.00 29.93           C  
ATOM   2828  O   ASP B 195      48.119  64.820  30.172  1.00 29.57           O  
ATOM   2829  CB  ASP B 195      51.229  65.163  30.017  1.00 30.08           C  
ATOM   2830  CG  ASP B 195      51.528  65.923  28.740  1.00 29.93           C  
ATOM   2831  OD1 ASP B 195      50.690  65.896  27.812  1.00 29.45           O  
ATOM   2832  OD2 ASP B 195      52.610  66.542  28.657  1.00 30.52           O  
ATOM   2833  N   SER B 196      48.410  66.988  29.641  1.00 29.77           N  
ATOM   2834  CA  SER B 196      47.083  67.169  29.025  1.00 29.65           C  
ATOM   2835  C   SER B 196      46.856  66.405  27.704  1.00 29.37           C  
ATOM   2836  O   SER B 196      45.712  66.189  27.289  1.00 29.28           O  
ATOM   2837  CB  SER B 196      46.787  68.667  28.830  1.00 29.76           C  
ATOM   2838  OG  SER B 196      47.679  69.266  27.897  1.00 29.98           O  
ATOM   2839  N   HIS B 197      47.945  66.006  27.052  1.00 29.07           N  
ATOM   2840  CA  HIS B 197      47.869  65.265  25.793  1.00 28.81           C  
ATOM   2841  C   HIS B 197      47.810  63.753  26.007  1.00 28.51           C  
ATOM   2842  O   HIS B 197      47.579  62.993  25.063  1.00 28.53           O  
ATOM   2843  CB  HIS B 197      49.043  65.640  24.886  1.00 28.84           C  
ATOM   2844  CG  HIS B 197      49.211  67.116  24.711  1.00 29.11           C  
ATOM   2845  ND1 HIS B 197      50.029  67.872  25.522  1.00 29.10           N  
ATOM   2846  CD2 HIS B 197      48.644  67.980  23.836  1.00 29.17           C  
ATOM   2847  CE1 HIS B 197      49.972  69.137  25.144  1.00 29.35           C  
ATOM   2848  NE2 HIS B 197      49.138  69.229  24.124  1.00 29.44           N  
ATOM   2849  N   LYS B 198      48.002  63.324  27.252  1.00 28.03           N  
ATOM   2850  CA  LYS B 198      47.969  61.907  27.593  1.00 27.64           C  
ATOM   2851  C   LYS B 198      46.657  61.521  28.273  1.00 27.39           C  
ATOM   2852  O   LYS B 198      45.939  62.377  28.791  1.00 27.40           O  
ATOM   2853  CB  LYS B 198      49.151  61.542  28.496  1.00 27.60           C  
ATOM   2854  CG  LYS B 198      50.519  61.830  27.901  1.00 27.23           C  
ATOM   2855  CD  LYS B 198      51.549  60.879  28.491  1.00 27.03           C  
ATOM   2856  CE  LYS B 198      52.921  61.045  27.865  1.00 25.74           C  
ATOM   2857  NZ  LYS B 198      53.657  62.171  28.476  1.00 25.70           N  
ATOM   2858  N   ASN B 199      46.347  60.228  28.251  1.00 27.04           N  
ATOM   2859  CA  ASN B 199      45.213  59.691  28.992  1.00 26.76           C  
ATOM   2860  C   ASN B 199      45.715  59.006  30.243  1.00 26.83           C  
ATOM   2861  O   ASN B 199      46.781  58.383  30.242  1.00 27.15           O  
ATOM   2862  CB  ASN B 199      44.430  58.675  28.165  1.00 26.66           C  
ATOM   2863  CG  ASN B 199      44.197  59.128  26.752  1.00 26.23           C  
ATOM   2864  OD1 ASN B 199      43.537  60.135  26.509  1.00 25.95           O  
ATOM   2865  ND2 ASN B 199      44.735  58.379  25.803  1.00 25.45           N  
ATOM   2866  N   TYR B 200      44.950  59.127  31.313  1.00 26.63           N  
ATOM   2867  CA  TYR B 200      45.283  58.463  32.547  1.00 26.69           C  
ATOM   2868  C   TYR B 200      44.032  57.802  33.090  1.00 26.89           C  
ATOM   2869  O   TYR B 200      42.915  58.252  32.825  1.00 26.74           O  
ATOM   2870  CB  TYR B 200      45.832  59.462  33.564  1.00 26.71           C  
ATOM   2871  CG  TYR B 200      47.178  60.065  33.218  1.00 26.22           C  
ATOM   2872  CD1 TYR B 200      48.360  59.520  33.723  1.00 26.43           C  
ATOM   2873  CD2 TYR B 200      47.271  61.198  32.411  1.00 26.65           C  
ATOM   2874  CE1 TYR B 200      49.606  60.082  33.425  1.00 26.00           C  
ATOM   2875  CE2 TYR B 200      48.511  61.767  32.100  1.00 26.66           C  
ATOM   2876  CZ  TYR B 200      49.672  61.204  32.611  1.00 26.57           C  
ATOM   2877  OH  TYR B 200      50.891  61.763  32.305  1.00 26.50           O  
ATOM   2878  N   VAL B 201      44.228  56.709  33.817  1.00 27.22           N  
ATOM   2879  CA  VAL B 201      43.169  56.129  34.629  1.00 27.65           C  
ATOM   2880  C   VAL B 201      43.643  56.076  36.076  1.00 27.51           C  
ATOM   2881  O   VAL B 201      44.820  55.859  36.350  1.00 27.60           O  
ATOM   2882  CB  VAL B 201      42.705  54.739  34.101  1.00 28.02           C  
ATOM   2883  CG1 VAL B 201      43.856  53.721  34.085  1.00 28.31           C  
ATOM   2884  CG2 VAL B 201      41.506  54.221  34.893  1.00 28.39           C  
ATOM   2885  N   VAL B 202      42.724  56.322  36.997  1.00 27.65           N  
ATOM   2886  CA  VAL B 202      43.052  56.329  38.416  1.00 27.20           C  
ATOM   2887  C   VAL B 202      42.205  55.313  39.123  1.00 27.07           C  
ATOM   2888  O   VAL B 202      40.976  55.324  39.016  1.00 27.05           O  
ATOM   2889  CB  VAL B 202      42.799  57.692  39.065  1.00 27.19           C  
ATOM   2890  CG1 VAL B 202      42.987  57.597  40.587  1.00 27.56           C  
ATOM   2891  CG2 VAL B 202      43.716  58.753  38.458  1.00 26.56           C  
ATOM   2892  N   ILE B 203      42.876  54.417  39.832  1.00 26.96           N  
ATOM   2893  CA  ILE B 203      42.189  53.490  40.702  1.00 26.70           C  
ATOM   2894  C   ILE B 203      42.139  54.106  42.097  1.00 27.02           C  
ATOM   2895  O   ILE B 203      43.133  54.658  42.585  1.00 27.04           O  
ATOM   2896  CB  ILE B 203      42.812  52.064  40.659  1.00 26.56           C  
ATOM   2897  CG1 ILE B 203      41.986  51.082  41.510  1.00 25.79           C  
ATOM   2898  CG2 ILE B 203      44.317  52.080  41.008  1.00 26.28           C  
ATOM   2899  CD1 ILE B 203      41.966  49.673  40.973  1.00 23.32           C  
ATOM   2900  N   VAL B 204      40.959  54.045  42.706  1.00 27.10           N  
ATOM   2901  CA  VAL B 204      40.709  54.667  43.991  1.00 27.37           C  
ATOM   2902  C   VAL B 204      40.298  53.585  44.984  1.00 27.74           C  
ATOM   2903  O   VAL B 204      39.375  52.814  44.726  1.00 28.14           O  
ATOM   2904  CB  VAL B 204      39.593  55.745  43.874  1.00 27.48           C  
ATOM   2905  CG1 VAL B 204      39.355  56.438  45.201  1.00 27.18           C  
ATOM   2906  CG2 VAL B 204      39.940  56.773  42.795  1.00 27.16           C  
ATOM   2907  N   ASN B 205      41.002  53.505  46.107  1.00 28.07           N  
ATOM   2908  CA  ASN B 205      40.549  52.683  47.225  1.00 28.01           C  
ATOM   2909  C   ASN B 205      40.010  53.609  48.309  1.00 28.02           C  
ATOM   2910  O   ASN B 205      40.771  54.155  49.110  1.00 27.78           O  
ATOM   2911  CB  ASN B 205      41.668  51.783  47.754  1.00 27.91           C  
ATOM   2912  CG  ASN B 205      41.177  50.770  48.793  1.00 28.57           C  
ATOM   2913  OD1 ASN B 205      40.045  50.852  49.286  1.00 29.58           O  
ATOM   2914  ND2 ASN B 205      42.036  49.813  49.136  1.00 27.37           N  
ATOM   2915  N   ASP B 206      38.690  53.794  48.303  1.00 28.15           N  
ATOM   2916  CA  ASP B 206      38.031  54.674  49.250  1.00 28.36           C  
ATOM   2917  C   ASP B 206      37.705  53.933  50.541  1.00 28.47           C  
ATOM   2918  O   ASP B 206      36.700  53.217  50.627  1.00 28.93           O  
ATOM   2919  CB  ASP B 206      36.763  55.263  48.635  1.00 28.54           C  
ATOM   2920  CG  ASP B 206      36.256  56.481  49.392  1.00 29.58           C  
ATOM   2921  OD1 ASP B 206      36.155  56.437  50.644  1.00 30.81           O  
ATOM   2922  OD2 ASP B 206      35.949  57.496  48.730  1.00 30.25           O  
ATOM   2923  N   GLY B 207      38.558  54.117  51.544  1.00 28.16           N  
ATOM   2924  CA  GLY B 207      38.389  53.470  52.831  1.00 27.68           C  
ATOM   2925  C   GLY B 207      37.222  53.971  53.651  1.00 27.65           C  
ATOM   2926  O   GLY B 207      36.803  53.302  54.594  1.00 27.92           O  
ATOM   2927  N   ARG B 208      36.696  55.145  53.313  1.00 27.42           N  
ATOM   2928  CA  ARG B 208      35.510  55.670  53.994  1.00 27.40           C  
ATOM   2929  C   ARG B 208      34.244  54.981  53.498  1.00 27.26           C  
ATOM   2930  O   ARG B 208      33.377  54.614  54.294  1.00 27.12           O  
ATOM   2931  CB  ARG B 208      35.384  57.191  53.826  1.00 27.61           C  
ATOM   2932  CG  ARG B 208      36.532  57.980  54.416  1.00 27.85           C  
ATOM   2933  CD  ARG B 208      37.668  58.112  53.416  1.00 29.92           C  
ATOM   2934  NE  ARG B 208      37.790  59.471  52.899  1.00 32.23           N  
ATOM   2935  CZ  ARG B 208      37.354  59.882  51.710  1.00 33.45           C  
ATOM   2936  NH1 ARG B 208      36.749  59.044  50.878  1.00 33.88           N  
ATOM   2937  NH2 ARG B 208      37.525  61.148  51.355  1.00 33.84           N  
ATOM   2938  N   LEU B 209      34.161  54.801  52.181  1.00 27.04           N  
ATOM   2939  CA  LEU B 209      33.023  54.140  51.551  1.00 26.92           C  
ATOM   2940  C   LEU B 209      33.148  52.639  51.608  1.00 26.64           C  
ATOM   2941  O   LEU B 209      32.140  51.944  51.494  1.00 26.87           O  
ATOM   2942  CB  LEU B 209      32.912  54.537  50.083  1.00 26.96           C  
ATOM   2943  CG  LEU B 209      32.305  55.850  49.610  1.00 26.82           C  
ATOM   2944  CD1 LEU B 209      32.918  57.019  50.317  1.00 29.10           C  
ATOM   2945  CD2 LEU B 209      32.575  55.975  48.129  1.00 27.26           C  
ATOM   2946  N   GLY B 210      34.383  52.151  51.770  1.00 26.42           N  
ATOM   2947  CA  GLY B 210      34.713  50.719  51.647  1.00 25.75           C  
ATOM   2948  C   GLY B 210      34.512  50.255  50.217  1.00 25.83           C  
ATOM   2949  O   GLY B 210      33.916  49.208  49.972  1.00 25.78           O  
ATOM   2950  N   HIS B 211      35.023  51.046  49.272  1.00 25.67           N  
ATOM   2951  CA  HIS B 211      34.646  50.936  47.867  1.00 25.47           C  
ATOM   2952  C   HIS B 211      35.809  51.224  46.910  1.00 25.70           C  
ATOM   2953  O   HIS B 211      36.490  52.253  47.026  1.00 25.78           O  
ATOM   2954  CB  HIS B 211      33.474  51.885  47.576  1.00 24.88           C  
ATOM   2955  CG  HIS B 211      32.922  51.772  46.187  1.00 24.28           C  
ATOM   2956  ND1 HIS B 211      32.416  50.595  45.678  1.00 22.35           N  
ATOM   2957  CD2 HIS B 211      32.777  52.699  45.208  1.00 23.52           C  
ATOM   2958  CE1 HIS B 211      31.994  50.800  44.443  1.00 22.72           C  
ATOM   2959  NE2 HIS B 211      32.196  52.069  44.135  1.00 22.05           N  
ATOM   2960  N   LYS B 212      36.031  50.304  45.977  1.00 25.55           N  
ATOM   2961  CA  LYS B 212      36.938  50.556  44.863  1.00 26.02           C  
ATOM   2962  C   LYS B 212      36.153  51.068  43.658  1.00 25.70           C  
ATOM   2963  O   LYS B 212      35.003  50.667  43.448  1.00 25.76           O  
ATOM   2964  CB  LYS B 212      37.752  49.307  44.493  1.00 26.04           C  
ATOM   2965  CG  LYS B 212      38.875  48.997  45.465  1.00 26.87           C  
ATOM   2966  CD  LYS B 212      39.900  48.068  44.860  1.00 29.49           C  
ATOM   2967  CE  LYS B 212      41.022  47.741  45.855  1.00 31.72           C  
ATOM   2968  NZ  LYS B 212      40.641  46.688  46.844  1.00 33.30           N  
ATOM   2969  N   PHE B 213      36.773  51.974  42.900  1.00 25.28           N  
ATOM   2970  CA  PHE B 213      36.232  52.459  41.632  1.00 25.17           C  
ATOM   2971  C   PHE B 213      37.330  53.137  40.818  1.00 25.27           C  
ATOM   2972  O   PHE B 213      38.467  53.277  41.291  1.00 25.04           O  
ATOM   2973  CB  PHE B 213      35.007  53.373  41.835  1.00 25.42           C  
ATOM   2974  CG  PHE B 213      35.310  54.696  42.515  1.00 25.66           C  
ATOM   2975  CD1 PHE B 213      35.423  55.870  41.764  1.00 25.39           C  
ATOM   2976  CD2 PHE B 213      35.459  54.770  43.897  1.00 24.80           C  
ATOM   2977  CE1 PHE B 213      35.691  57.102  42.379  1.00 25.50           C  
ATOM   2978  CE2 PHE B 213      35.731  55.992  44.521  1.00 26.11           C  
ATOM   2979  CZ  PHE B 213      35.840  57.165  43.758  1.00 25.59           C  
ATOM   2980  N   LEU B 214      36.986  53.550  39.595  1.00 25.42           N  
ATOM   2981  CA  LEU B 214      37.964  54.111  38.659  1.00 25.13           C  
ATOM   2982  C   LEU B 214      37.570  55.488  38.152  1.00 25.03           C  
ATOM   2983  O   LEU B 214      36.394  55.750  37.876  1.00 25.08           O  
ATOM   2984  CB  LEU B 214      38.175  53.179  37.459  1.00 25.08           C  
ATOM   2985  CG  LEU B 214      38.575  51.712  37.668  1.00 25.13           C  
ATOM   2986  CD1 LEU B 214      38.536  50.991  36.351  1.00 24.79           C  
ATOM   2987  CD2 LEU B 214      39.948  51.575  38.274  1.00 25.36           C  
ATOM   2988  N   ILE B 215      38.569  56.362  38.037  1.00 24.88           N  
ATOM   2989  CA  ILE B 215      38.398  57.658  37.389  1.00 24.62           C  
ATOM   2990  C   ILE B 215      39.172  57.673  36.076  1.00 24.68           C  
ATOM   2991  O   ILE B 215      40.396  57.496  36.062  1.00 24.46           O  
ATOM   2992  CB  ILE B 215      38.853  58.842  38.282  1.00 24.74           C  
ATOM   2993  CG1 ILE B 215      38.183  58.774  39.665  1.00 24.18           C  
ATOM   2994  CG2 ILE B 215      38.554  60.169  37.589  1.00 23.32           C  
ATOM   2995  CD1 ILE B 215      38.745  59.755  40.680  1.00 24.45           C  
ATOM   2996  N   ASP B 216      38.440  57.871  34.982  1.00 24.70           N  
ATOM   2997  CA  ASP B 216      39.018  57.881  33.644  1.00 24.97           C  
ATOM   2998  C   ASP B 216      39.333  59.318  33.231  1.00 25.02           C  
ATOM   2999  O   ASP B 216      38.461  60.199  33.295  1.00 24.92           O  
ATOM   3000  CB  ASP B 216      38.065  57.219  32.642  1.00 25.11           C  
ATOM   3001  CG  ASP B 216      38.741  56.878  31.314  1.00 26.13           C  
ATOM   3002  OD1 ASP B 216      39.963  57.109  31.171  1.00 27.15           O  
ATOM   3003  OD2 ASP B 216      38.049  56.365  30.404  1.00 27.38           O  
ATOM   3004  N   LEU B 217      40.581  59.550  32.820  1.00 24.89           N  
ATOM   3005  CA  LEU B 217      41.026  60.896  32.465  1.00 25.30           C  
ATOM   3006  C   LEU B 217      41.549  60.968  31.031  1.00 25.40           C  
ATOM   3007  O   LEU B 217      42.754  61.017  30.814  1.00 25.57           O  
ATOM   3008  CB  LEU B 217      42.077  61.412  33.465  1.00 25.01           C  
ATOM   3009  CG  LEU B 217      41.721  61.441  34.960  1.00 25.23           C  
ATOM   3010  CD1 LEU B 217      42.931  61.852  35.780  1.00 26.09           C  
ATOM   3011  CD2 LEU B 217      40.537  62.351  35.275  1.00 25.47           C  
ATOM   3012  N   PRO B 218      40.641  60.959  30.040  1.00 25.71           N  
ATOM   3013  CA  PRO B 218      41.079  61.052  28.646  1.00 26.16           C  
ATOM   3014  C   PRO B 218      41.691  62.407  28.304  1.00 26.56           C  
ATOM   3015  O   PRO B 218      41.351  63.418  28.915  1.00 27.08           O  
ATOM   3016  CB  PRO B 218      39.780  60.853  27.860  1.00 26.10           C  
ATOM   3017  CG  PRO B 218      38.833  60.202  28.837  1.00 26.11           C  
ATOM   3018  CD  PRO B 218      39.180  60.822  30.137  1.00 25.68           C  
ATOM   3019  N   ALA B 219      42.600  62.420  27.338  1.00 26.81           N  
ATOM   3020  CA  ALA B 219      43.150  63.666  26.820  1.00 26.71           C  
ATOM   3021  C   ALA B 219      42.086  64.394  25.995  1.00 26.74           C  
ATOM   3022  O   ALA B 219      41.366  63.775  25.199  1.00 26.41           O  
ATOM   3023  CB  ALA B 219      44.381  63.380  25.979  1.00 26.70           C  
ATOM   3024  N   PRO B 224      40.701  72.193  30.233  1.00 34.78           N  
ATOM   3025  CA  PRO B 224      39.368  71.812  30.697  1.00 34.64           C  
ATOM   3026  C   PRO B 224      39.214  70.293  30.780  1.00 34.65           C  
ATOM   3027  O   PRO B 224      38.353  69.714  30.102  1.00 34.86           O  
ATOM   3028  CB  PRO B 224      38.432  72.402  29.627  1.00 34.62           C  
ATOM   3029  CG  PRO B 224      39.341  72.956  28.517  1.00 34.82           C  
ATOM   3030  CD  PRO B 224      40.730  72.475  28.788  1.00 34.72           C  
ATOM   3031  N   ARG B 225      40.049  69.664  31.609  1.00 34.39           N  
ATOM   3032  CA  ARG B 225      40.062  68.204  31.774  1.00 34.14           C  
ATOM   3033  C   ARG B 225      38.722  67.635  32.230  1.00 33.69           C  
ATOM   3034  O   ARG B 225      38.026  68.244  33.039  1.00 33.83           O  
ATOM   3035  CB  ARG B 225      41.138  67.787  32.769  1.00 34.29           C  
ATOM   3036  CG  ARG B 225      42.441  67.350  32.149  1.00 34.52           C  
ATOM   3037  CD  ARG B 225      43.089  66.339  33.067  1.00 34.84           C  
ATOM   3038  NE  ARG B 225      44.445  65.992  32.664  1.00 35.24           N  
ATOM   3039  CZ  ARG B 225      44.757  65.006  31.828  1.00 35.02           C  
ATOM   3040  NH1 ARG B 225      43.814  64.251  31.270  1.00 33.63           N  
ATOM   3041  NH2 ARG B 225      46.029  64.780  31.549  1.00 35.88           N  
ATOM   3042  N   THR B 226      38.376  66.460  31.716  1.00 33.07           N  
ATOM   3043  CA  THR B 226      37.073  65.853  31.986  1.00 32.59           C  
ATOM   3044  C   THR B 226      37.253  64.475  32.611  1.00 31.96           C  
ATOM   3045  O   THR B 226      38.153  63.727  32.230  1.00 32.22           O  
ATOM   3046  CB  THR B 226      36.237  65.723  30.691  1.00 32.79           C  
ATOM   3047  OG1 THR B 226      36.966  64.947  29.734  1.00 33.12           O  
ATOM   3048  CG2 THR B 226      35.923  67.110  30.082  1.00 32.97           C  
ATOM   3049  N   ALA B 227      36.393  64.133  33.562  1.00 31.05           N  
ATOM   3050  CA  ALA B 227      36.541  62.885  34.299  1.00 30.06           C  
ATOM   3051  C   ALA B 227      35.302  62.007  34.192  1.00 29.72           C  
ATOM   3052  O   ALA B 227      34.177  62.506  34.160  1.00 29.51           O  
ATOM   3053  CB  ALA B 227      36.871  63.172  35.743  1.00 29.85           C  
ATOM   3054  N   TYR B 228      35.520  60.693  34.126  1.00 29.16           N  
ATOM   3055  CA  TYR B 228      34.430  59.728  34.075  1.00 28.45           C  
ATOM   3056  C   TYR B 228      34.644  58.674  35.143  1.00 28.40           C  
ATOM   3057  O   TYR B 228      35.783  58.325  35.471  1.00 28.55           O  
ATOM   3058  CB  TYR B 228      34.359  59.047  32.704  1.00 28.22           C  
ATOM   3059  CG  TYR B 228      34.165  59.979  31.530  1.00 27.61           C  
ATOM   3060  CD1 TYR B 228      32.898  60.441  31.185  1.00 26.42           C  
ATOM   3061  CD2 TYR B 228      35.253  60.384  30.751  1.00 27.56           C  
ATOM   3062  CE1 TYR B 228      32.713  61.289  30.105  1.00 26.01           C  
ATOM   3063  CE2 TYR B 228      35.080  61.234  29.664  1.00 26.79           C  
ATOM   3064  CZ  TYR B 228      33.805  61.680  29.351  1.00 27.06           C  
ATOM   3065  OH  TYR B 228      33.619  62.521  28.283  1.00 27.57           O  
ATOM   3066  N   ILE B 229      33.545  58.163  35.683  1.00 28.23           N  
ATOM   3067  CA  ILE B 229      33.612  57.090  36.671  1.00 28.08           C  
ATOM   3068  C   ILE B 229      33.257  55.742  36.043  1.00 27.88           C  
ATOM   3069  O   ILE B 229      32.341  55.648  35.217  1.00 27.82           O  
ATOM   3070  CB  ILE B 229      32.705  57.379  37.904  1.00 28.04           C  
ATOM   3071  CG1 ILE B 229      33.174  58.642  38.658  1.00 28.36           C  
ATOM   3072  CG2 ILE B 229      32.624  56.167  38.820  1.00 27.36           C  
ATOM   3073  CD1 ILE B 229      34.662  58.683  39.055  1.00 26.76           C  
ATOM   3074  N   ILE B 230      34.031  54.726  36.416  1.00 27.67           N  
ATOM   3075  CA  ILE B 230      33.765  53.334  36.092  1.00 27.44           C  
ATOM   3076  C   ILE B 230      33.784  52.586  37.415  1.00 27.55           C  
ATOM   3077  O   ILE B 230      34.784  52.618  38.138  1.00 27.58           O  
ATOM   3078  CB  ILE B 230      34.844  52.728  35.156  1.00 27.48           C  
ATOM   3079  CG1 ILE B 230      34.891  53.466  33.809  1.00 27.50           C  
ATOM   3080  CG2 ILE B 230      34.580  51.242  34.932  1.00 27.25           C  
ATOM   3081  CD1 ILE B 230      36.038  53.012  32.879  1.00 27.14           C  
ATOM   3082  N   GLN B 231      32.677  51.927  37.735  1.00 27.85           N  
ATOM   3083  CA  GLN B 231      32.549  51.193  38.986  1.00 28.34           C  
ATOM   3084  C   GLN B 231      31.508  50.086  38.892  1.00 28.76           C  
ATOM   3085  O   GLN B 231      30.692  50.065  37.970  1.00 28.71           O  
ATOM   3086  CB  GLN B 231      32.137  52.152  40.108  1.00 28.36           C  
ATOM   3087  CG  GLN B 231      30.762  52.792  39.905  1.00 28.98           C  
ATOM   3088  CD  GLN B 231      30.219  53.448  41.158  1.00 28.61           C  
ATOM   3089  OE1 GLN B 231      30.721  53.222  42.253  1.00 29.91           O  
ATOM   3090  NE2 GLN B 231      29.184  54.262  41.000  1.00 27.87           N  
ATOM   3091  N   SER B 232      31.549  49.163  39.848  1.00 29.27           N  
ATOM   3092  CA  SER B 232      30.375  48.366  40.183  1.00 30.07           C  
ATOM   3093  C   SER B 232      30.278  48.323  41.702  1.00 30.55           C  
ATOM   3094  O   SER B 232      31.227  48.689  42.387  1.00 30.94           O  
ATOM   3095  CB  SER B 232      30.457  46.969  39.595  1.00 29.84           C  
ATOM   3096  OG  SER B 232      31.434  46.227  40.272  1.00 30.91           O  
ATOM   3097  N   ASP B 233      29.128  47.920  42.228  1.00 30.81           N  
ATOM   3098  CA  ASP B 233      28.931  47.866  43.662  1.00 31.22           C  
ATOM   3099  C   ASP B 233      27.831  46.850  43.904  1.00 31.32           C  
ATOM   3100  O   ASP B 233      26.802  46.893  43.225  1.00 31.64           O  
ATOM   3101  CB  ASP B 233      28.516  49.247  44.204  1.00 31.53           C  
ATOM   3102  CG  ASP B 233      28.782  49.405  45.714  1.00 32.49           C  
ATOM   3103  OD1 ASP B 233      28.570  48.436  46.486  1.00 32.77           O  
ATOM   3104  OD2 ASP B 233      29.202  50.508  46.132  1.00 33.54           O  
ATOM   3105  N   LEU B 234      28.037  45.936  44.854  1.00 31.02           N  
ATOM   3106  CA  LEU B 234      26.993  44.968  45.180  1.00 30.69           C  
ATOM   3107  C   LEU B 234      26.020  45.513  46.213  1.00 30.75           C  
ATOM   3108  O   LEU B 234      25.125  44.796  46.657  1.00 30.67           O  
ATOM   3109  CB  LEU B 234      27.573  43.630  45.640  1.00 30.64           C  
ATOM   3110  CG  LEU B 234      28.085  42.597  44.620  1.00 30.76           C  
ATOM   3111  CD1 LEU B 234      28.097  41.219  45.265  1.00 30.11           C  
ATOM   3112  CD2 LEU B 234      27.305  42.559  43.289  1.00 30.35           C  
ATOM   3113  N   GLY B 235      26.189  46.789  46.570  1.00 30.93           N  
ATOM   3114  CA  GLY B 235      25.369  47.453  47.583  1.00 30.75           C  
ATOM   3115  C   GLY B 235      25.689  46.973  48.986  1.00 31.04           C  
ATOM   3116  O   GLY B 235      26.446  46.016  49.169  1.00 30.78           O  
ATOM   3117  N   GLY B 236      25.105  47.640  49.979  1.00 31.31           N  
ATOM   3118  CA  GLY B 236      25.292  47.269  51.375  1.00 31.55           C  
ATOM   3119  C   GLY B 236      25.913  48.368  52.208  1.00 31.79           C  
ATOM   3120  O   GLY B 236      25.680  48.446  53.414  1.00 31.87           O  
ATOM   3121  N   GLY B 237      26.722  49.202  51.568  1.00 31.97           N  
ATOM   3122  CA  GLY B 237      27.342  50.344  52.232  1.00 32.41           C  
ATOM   3123  C   GLY B 237      26.605  51.616  51.874  1.00 32.67           C  
ATOM   3124  O   GLY B 237      25.405  51.580  51.567  1.00 33.12           O  
ATOM   3125  N   ALA B 238      27.318  52.739  51.886  1.00 32.59           N  
ATOM   3126  CA  ALA B 238      26.694  54.040  51.627  1.00 32.38           C  
ATOM   3127  C   ALA B 238      25.961  54.073  50.290  1.00 32.45           C  
ATOM   3128  O   ALA B 238      24.956  54.772  50.140  1.00 32.62           O  
ATOM   3129  CB  ALA B 238      27.724  55.156  51.697  1.00 32.02           C  
ATOM   3130  N   LEU B 239      26.455  53.301  49.330  1.00 32.46           N  
ATOM   3131  CA  LEU B 239      25.946  53.364  47.966  1.00 32.52           C  
ATOM   3132  C   LEU B 239      25.003  52.205  47.636  1.00 32.21           C  
ATOM   3133  O   LEU B 239      25.148  51.116  48.191  1.00 32.10           O  
ATOM   3134  CB  LEU B 239      27.114  53.431  46.972  1.00 32.81           C  
ATOM   3135  CG  LEU B 239      28.060  54.644  47.108  1.00 33.54           C  
ATOM   3136  CD1 LEU B 239      29.365  54.469  46.310  1.00 33.70           C  
ATOM   3137  CD2 LEU B 239      27.366  55.945  46.723  1.00 32.97           C  
ATOM   3138  N   PRO B 240      24.012  52.451  46.754  1.00 31.87           N  
ATOM   3139  CA  PRO B 240      23.165  51.393  46.189  1.00 31.46           C  
ATOM   3140  C   PRO B 240      23.934  50.482  45.238  1.00 30.93           C  
ATOM   3141  O   PRO B 240      24.882  50.925  44.580  1.00 30.83           O  
ATOM   3142  CB  PRO B 240      22.096  52.154  45.380  1.00 31.22           C  
ATOM   3143  CG  PRO B 240      22.358  53.551  45.519  1.00 31.53           C  
ATOM   3144  CD  PRO B 240      23.642  53.786  46.255  1.00 32.02           C  
ATOM   3145  N   ALA B 241      23.505  49.226  45.162  1.00 30.40           N  
ATOM   3146  CA  ALA B 241      24.061  48.264  44.223  1.00 30.06           C  
ATOM   3147  C   ALA B 241      23.868  48.750  42.794  1.00 29.99           C  
ATOM   3148  O   ALA B 241      22.823  49.312  42.453  1.00 29.90           O  
ATOM   3149  CB  ALA B 241      23.404  46.916  44.405  1.00 29.78           C  
ATOM   3150  N   VAL B 242      24.895  48.544  41.976  1.00 29.76           N  
ATOM   3151  CA  VAL B 242      24.845  48.865  40.563  1.00 29.48           C  
ATOM   3152  C   VAL B 242      25.753  47.913  39.804  1.00 29.87           C  
ATOM   3153  O   VAL B 242      26.913  47.694  40.179  1.00 29.63           O  
ATOM   3154  CB  VAL B 242      25.188  50.353  40.282  1.00 29.42           C  
ATOM   3155  CG1 VAL B 242      26.658  50.665  40.558  1.00 29.14           C  
ATOM   3156  CG2 VAL B 242      24.800  50.735  38.868  1.00 28.53           C  
ATOM   3157  N   ARG B 243      25.190  47.323  38.755  1.00 30.32           N  
ATOM   3158  CA  ARG B 243      25.920  46.420  37.886  1.00 30.68           C  
ATOM   3159  C   ARG B 243      26.784  47.254  36.952  1.00 30.86           C  
ATOM   3160  O   ARG B 243      26.321  48.263  36.404  1.00 30.88           O  
ATOM   3161  CB  ARG B 243      24.941  45.575  37.082  1.00 30.80           C  
ATOM   3162  CG  ARG B 243      24.070  44.650  37.921  1.00 32.01           C  
ATOM   3163  CD  ARG B 243      22.926  44.053  37.108  1.00 34.58           C  
ATOM   3164  NE  ARG B 243      23.405  43.189  36.031  1.00 37.16           N  
ATOM   3165  CZ  ARG B 243      23.462  43.541  34.746  1.00 39.63           C  
ATOM   3166  NH1 ARG B 243      23.061  44.751  34.357  1.00 40.65           N  
ATOM   3167  NH2 ARG B 243      23.923  42.681  33.842  1.00 39.72           N  
ATOM   3168  N   VAL B 244      28.035  46.839  36.779  1.00 30.83           N  
ATOM   3169  CA  VAL B 244      28.977  47.550  35.912  1.00 31.06           C  
ATOM   3170  C   VAL B 244      28.366  47.924  34.547  1.00 31.38           C  
ATOM   3171  O   VAL B 244      28.613  49.018  34.028  1.00 31.46           O  
ATOM   3172  CB  VAL B 244      30.343  46.767  35.748  1.00 31.18           C  
ATOM   3173  CG1 VAL B 244      30.150  45.408  35.052  1.00 30.54           C  
ATOM   3174  CG2 VAL B 244      31.379  47.605  35.012  1.00 30.72           C  
ATOM   3175  N   GLU B 245      27.561  47.019  33.985  1.00 31.85           N  
ATOM   3176  CA  GLU B 245      26.970  47.214  32.651  1.00 31.97           C  
ATOM   3177  C   GLU B 245      25.984  48.365  32.633  1.00 31.62           C  
ATOM   3178  O   GLU B 245      25.941  49.120  31.666  1.00 31.90           O  
ATOM   3179  CB  GLU B 245      26.277  45.948  32.144  1.00 32.02           C  
ATOM   3180  CG  GLU B 245      27.222  44.829  31.767  1.00 33.72           C  
ATOM   3181  CD  GLU B 245      27.322  43.743  32.832  1.00 36.98           C  
ATOM   3182  OE1 GLU B 245      26.909  43.986  33.996  1.00 37.46           O  
ATOM   3183  OE2 GLU B 245      27.818  42.636  32.493  1.00 37.63           O  
ATOM   3184  N   ASP B 246      25.200  48.487  33.701  1.00 31.26           N  
ATOM   3185  CA  ASP B 246      24.208  49.551  33.822  1.00 31.08           C  
ATOM   3186  C   ASP B 246      24.875  50.864  34.171  1.00 30.83           C  
ATOM   3187  O   ASP B 246      24.449  51.919  33.709  1.00 30.91           O  
ATOM   3188  CB  ASP B 246      23.166  49.210  34.884  1.00 31.13           C  
ATOM   3189  CG  ASP B 246      22.442  47.917  34.591  1.00 31.47           C  
ATOM   3190  OD1 ASP B 246      22.137  47.644  33.408  1.00 31.16           O  
ATOM   3191  OD2 ASP B 246      22.183  47.171  35.554  1.00 32.18           O  
ATOM   3192  N   TRP B 247      25.918  50.802  34.990  1.00 30.63           N  
ATOM   3193  CA  TRP B 247      26.679  51.998  35.290  1.00 30.36           C  
ATOM   3194  C   TRP B 247      27.273  52.592  34.022  1.00 30.36           C  
ATOM   3195  O   TRP B 247      27.084  53.778  33.751  1.00 30.36           O  
ATOM   3196  CB  TRP B 247      27.792  51.767  36.316  1.00 29.97           C  
ATOM   3197  CG  TRP B 247      28.533  53.058  36.517  1.00 30.02           C  
ATOM   3198  CD1 TRP B 247      29.678  53.461  35.882  1.00 29.15           C  
ATOM   3199  CD2 TRP B 247      28.123  54.154  37.343  1.00 29.58           C  
ATOM   3200  NE1 TRP B 247      30.014  54.727  36.286  1.00 29.08           N  
ATOM   3201  CE2 TRP B 247      29.083  55.176  37.185  1.00 28.97           C  
ATOM   3202  CE3 TRP B 247      27.039  54.367  38.215  1.00 29.24           C  
ATOM   3203  CZ2 TRP B 247      29.002  56.394  37.869  1.00 29.22           C  
ATOM   3204  CZ3 TRP B 247      26.957  55.579  38.897  1.00 29.16           C  
ATOM   3205  CH2 TRP B 247      27.934  56.577  38.718  1.00 29.27           C  
ATOM   3206  N   ILE B 248      27.998  51.772  33.261  1.00 30.59           N  
ATOM   3207  CA  ILE B 248      28.635  52.224  32.022  1.00 30.77           C  
ATOM   3208  C   ILE B 248      27.593  52.705  31.024  1.00 30.80           C  
ATOM   3209  O   ILE B 248      27.740  53.762  30.424  1.00 30.93           O  
ATOM   3210  CB  ILE B 248      29.540  51.133  31.421  1.00 30.96           C  
ATOM   3211  CG1 ILE B 248      30.816  51.003  32.264  1.00 31.60           C  
ATOM   3212  CG2 ILE B 248      29.892  51.441  29.965  1.00 30.70           C  
ATOM   3213  CD1 ILE B 248      31.672  49.801  31.929  1.00 31.97           C  
ATOM   3214  N   SER B 249      26.522  51.938  30.886  1.00 31.34           N  
ATOM   3215  CA  SER B 249      25.434  52.274  29.980  1.00 31.82           C  
ATOM   3216  C   SER B 249      24.850  53.644  30.286  1.00 31.99           C  
ATOM   3217  O   SER B 249      24.713  54.460  29.386  1.00 32.34           O  
ATOM   3218  CB  SER B 249      24.330  51.218  30.055  1.00 31.89           C  
ATOM   3219  OG  SER B 249      23.379  51.397  29.025  1.00 32.50           O  
ATOM   3220  N   ARG B 250      24.521  53.894  31.552  1.00 32.14           N  
ATOM   3221  CA  ARG B 250      23.799  55.111  31.938  1.00 32.35           C  
ATOM   3222  C   ARG B 250      24.704  56.303  32.244  1.00 32.32           C  
ATOM   3223  O   ARG B 250      24.328  57.445  32.002  1.00 32.10           O  
ATOM   3224  CB  ARG B 250      22.881  54.846  33.140  1.00 32.44           C  
ATOM   3225  CG  ARG B 250      21.765  53.821  32.900  1.00 32.97           C  
ATOM   3226  CD  ARG B 250      20.476  54.472  32.396  1.00 34.03           C  
ATOM   3227  NE  ARG B 250      19.984  55.507  33.306  1.00 34.69           N  
ATOM   3228  CZ  ARG B 250      19.360  55.272  34.461  1.00 35.39           C  
ATOM   3229  NH1 ARG B 250      19.139  54.031  34.879  1.00 35.39           N  
ATOM   3230  NH2 ARG B 250      18.964  56.287  35.213  1.00 35.02           N  
ATOM   3231  N   ARG B 251      25.893  56.042  32.773  1.00 32.55           N  
ATOM   3232  CA  ARG B 251      26.707  57.111  33.337  1.00 32.85           C  
ATOM   3233  C   ARG B 251      28.176  57.069  32.937  1.00 32.71           C  
ATOM   3234  O   ARG B 251      28.950  57.948  33.331  1.00 32.94           O  
ATOM   3235  CB  ARG B 251      26.573  57.124  34.867  1.00 33.28           C  
ATOM   3236  CG  ARG B 251      25.221  57.663  35.376  1.00 34.91           C  
ATOM   3237  CD  ARG B 251      25.372  58.402  36.692  1.00 38.01           C  
ATOM   3238  NE  ARG B 251      26.285  59.541  36.580  1.00 41.28           N  
ATOM   3239  CZ  ARG B 251      26.738  60.258  37.609  1.00 43.14           C  
ATOM   3240  NH1 ARG B 251      26.371  59.975  38.855  1.00 43.17           N  
ATOM   3241  NH2 ARG B 251      27.567  61.270  37.389  1.00 45.37           N  
ATOM   3242  N   GLY B 252      28.552  56.072  32.139  1.00 32.31           N  
ATOM   3243  CA  GLY B 252      29.947  55.870  31.755  1.00 32.02           C  
ATOM   3244  C   GLY B 252      30.511  56.954  30.860  1.00 32.00           C  
ATOM   3245  O   GLY B 252      31.725  57.159  30.820  1.00 31.64           O  
ATOM   3246  N   SER B 253      29.627  57.640  30.139  1.00 32.06           N  
ATOM   3247  CA  SER B 253      30.024  58.728  29.249  1.00 32.41           C  
ATOM   3248  C   SER B 253      29.568  60.064  29.818  1.00 32.58           C  
ATOM   3249  O   SER B 253      29.464  61.066  29.100  1.00 32.29           O  
ATOM   3250  CB  SER B 253      29.452  58.515  27.848  1.00 32.38           C  
ATOM   3251  OG  SER B 253      29.904  57.283  27.322  1.00 33.00           O  
ATOM   3252  N   ASP B 254      29.319  60.064  31.127  1.00 32.84           N  
ATOM   3253  CA  ASP B 254      28.837  61.242  31.823  1.00 33.02           C  
ATOM   3254  C   ASP B 254      29.961  61.913  32.624  1.00 32.72           C  
ATOM   3255  O   ASP B 254      30.386  61.403  33.653  1.00 32.48           O  
ATOM   3256  CB  ASP B 254      27.652  60.866  32.713  1.00 33.21           C  
ATOM   3257  CG  ASP B 254      26.797  62.054  33.067  1.00 34.58           C  
ATOM   3258  OD1 ASP B 254      26.704  62.990  32.242  1.00 35.10           O  
ATOM   3259  OD2 ASP B 254      26.220  62.055  34.178  1.00 37.70           O  
ATOM   3260  N   PRO B 255      30.449  63.065  32.143  1.00 32.81           N  
ATOM   3261  CA  PRO B 255      31.582  63.725  32.788  1.00 33.02           C  
ATOM   3262  C   PRO B 255      31.247  64.205  34.196  1.00 33.32           C  
ATOM   3263  O   PRO B 255      30.144  64.717  34.442  1.00 33.30           O  
ATOM   3264  CB  PRO B 255      31.865  64.931  31.881  1.00 32.86           C  
ATOM   3265  CG  PRO B 255      31.142  64.658  30.612  1.00 32.85           C  
ATOM   3266  CD  PRO B 255      29.962  63.826  30.982  1.00 32.92           C  
ATOM   3267  N   VAL B 256      32.193  64.024  35.111  1.00 33.39           N  
ATOM   3268  CA  VAL B 256      32.039  64.523  36.468  1.00 33.46           C  
ATOM   3269  C   VAL B 256      32.971  65.701  36.710  1.00 33.43           C  
ATOM   3270  O   VAL B 256      34.106  65.738  36.209  1.00 33.40           O  
ATOM   3271  CB  VAL B 256      32.240  63.421  37.544  1.00 33.45           C  
ATOM   3272  CG1 VAL B 256      31.177  62.332  37.395  1.00 33.56           C  
ATOM   3273  CG2 VAL B 256      33.634  62.827  37.463  1.00 33.53           C  
ATOM   3274  N   SER B 257      32.457  66.667  37.468  1.00 33.29           N  
ATOM   3275  CA  SER B 257      33.190  67.864  37.842  1.00 33.09           C  
ATOM   3276  C   SER B 257      34.478  67.482  38.556  1.00 32.83           C  
ATOM   3277  O   SER B 257      34.460  66.672  39.479  1.00 32.97           O  
ATOM   3278  CB  SER B 257      32.307  68.742  38.737  1.00 33.10           C  
ATOM   3279  OG  SER B 257      33.070  69.662  39.496  1.00 33.72           O  
ATOM   3280  N   LEU B 258      35.590  68.052  38.107  1.00 32.71           N  
ATOM   3281  CA  LEU B 258      36.889  67.861  38.756  1.00 32.62           C  
ATOM   3282  C   LEU B 258      36.934  68.434  40.168  1.00 32.43           C  
ATOM   3283  O   LEU B 258      37.657  67.930  41.027  1.00 32.36           O  
ATOM   3284  CB  LEU B 258      37.997  68.507  37.932  1.00 32.59           C  
ATOM   3285  CG  LEU B 258      38.495  67.743  36.713  1.00 32.85           C  
ATOM   3286  CD1 LEU B 258      39.626  68.536  36.065  1.00 33.12           C  
ATOM   3287  CD2 LEU B 258      38.947  66.341  37.107  1.00 32.15           C  
ATOM   3288  N   ASP B 259      36.167  69.499  40.380  1.00 32.16           N  
ATOM   3289  CA  ASP B 259      36.043  70.161  41.673  1.00 31.98           C  
ATOM   3290  C   ASP B 259      35.362  69.232  42.684  1.00 31.56           C  
ATOM   3291  O   ASP B 259      35.862  69.040  43.802  1.00 31.34           O  
ATOM   3292  CB  ASP B 259      35.247  71.460  41.489  1.00 32.28           C  
ATOM   3293  CG  ASP B 259      35.196  72.316  42.739  1.00 33.39           C  
ATOM   3294  OD1 ASP B 259      36.098  72.217  43.605  1.00 34.93           O  
ATOM   3295  OD2 ASP B 259      34.242  73.116  42.843  1.00 35.02           O  
ATOM   3296  N   GLU B 260      34.234  68.647  42.274  1.00 30.95           N  
ATOM   3297  CA  GLU B 260      33.441  67.767  43.137  1.00 30.53           C  
ATOM   3298  C   GLU B 260      34.200  66.510  43.504  1.00 30.10           C  
ATOM   3299  O   GLU B 260      34.116  66.020  44.637  1.00 29.88           O  
ATOM   3300  CB  GLU B 260      32.130  67.391  42.453  1.00 30.54           C  
ATOM   3301  CG  GLU B 260      31.093  68.476  42.553  1.00 31.23           C  
ATOM   3302  CD  GLU B 260      29.810  68.138  41.838  1.00 32.19           C  
ATOM   3303  OE1 GLU B 260      29.304  67.006  41.997  1.00 30.77           O  
ATOM   3304  OE2 GLU B 260      29.303  69.028  41.120  1.00 33.87           O  
ATOM   3305  N   LEU B 261      34.927  65.996  42.516  1.00 29.54           N  
ATOM   3306  CA  LEU B 261      35.766  64.831  42.666  1.00 29.01           C  
ATOM   3307  C   LEU B 261      36.877  65.177  43.642  1.00 29.03           C  
ATOM   3308  O   LEU B 261      37.184  64.404  44.548  1.00 28.84           O  
ATOM   3309  CB  LEU B 261      36.330  64.448  41.302  1.00 28.87           C  
ATOM   3310  CG  LEU B 261      36.606  63.001  40.893  1.00 28.50           C  
ATOM   3311  CD1 LEU B 261      35.464  62.036  41.222  1.00 27.17           C  
ATOM   3312  CD2 LEU B 261      36.874  63.001  39.402  1.00 28.53           C  
ATOM   3313  N   ASN B 262      37.452  66.363  43.479  1.00 29.16           N  
ATOM   3314  CA  ASN B 262      38.449  66.845  44.418  1.00 29.64           C  
ATOM   3315  C   ASN B 262      37.911  66.980  45.835  1.00 29.86           C  
ATOM   3316  O   ASN B 262      38.639  66.734  46.798  1.00 30.05           O  
ATOM   3317  CB  ASN B 262      39.051  68.170  43.969  1.00 29.50           C  
ATOM   3318  CG  ASN B 262      40.199  68.598  44.846  1.00 30.28           C  
ATOM   3319  OD1 ASN B 262      40.383  69.788  45.110  1.00 31.53           O  
ATOM   3320  ND2 ASN B 262      40.970  67.624  45.334  1.00 31.02           N  
ATOM   3321  N   GLN B 263      36.646  67.370  45.955  1.00 30.12           N  
ATOM   3322  CA  GLN B 263      35.987  67.451  47.251  1.00 30.69           C  
ATOM   3323  C   GLN B 263      35.807  66.086  47.912  1.00 30.61           C  
ATOM   3324  O   GLN B 263      36.076  65.940  49.103  1.00 30.80           O  
ATOM   3325  CB  GLN B 263      34.629  68.140  47.130  1.00 30.84           C  
ATOM   3326  CG  GLN B 263      34.685  69.651  46.968  1.00 31.23           C  
ATOM   3327  CD  GLN B 263      33.331  70.231  46.583  1.00 31.61           C  
ATOM   3328  OE1 GLN B 263      32.327  70.002  47.270  1.00 33.04           O  
ATOM   3329  NE2 GLN B 263      33.295  70.981  45.475  1.00 31.40           N  
ATOM   3330  N   LEU B 264      35.339  65.104  47.141  1.00 30.62           N  
ATOM   3331  CA  LEU B 264      35.106  63.741  47.635  1.00 30.82           C  
ATOM   3332  C   LEU B 264      36.407  63.066  48.072  1.00 31.17           C  
ATOM   3333  O   LEU B 264      36.487  62.507  49.166  1.00 31.28           O  
ATOM   3334  CB  LEU B 264      34.469  62.894  46.532  1.00 30.82           C  
ATOM   3335  CG  LEU B 264      33.641  61.623  46.758  1.00 30.16           C  
ATOM   3336  CD1 LEU B 264      33.862  60.714  45.559  1.00 28.46           C  
ATOM   3337  CD2 LEU B 264      33.925  60.873  48.051  1.00 29.06           C  
ATOM   3338  N   LEU B 265      37.418  63.125  47.205  1.00 31.30           N  
ATOM   3339  CA  LEU B 265      38.703  62.478  47.438  1.00 31.49           C  
ATOM   3340  C   LEU B 265      39.589  63.239  48.419  1.00 31.68           C  
ATOM   3341  O   LEU B 265      40.786  62.964  48.529  1.00 32.15           O  
ATOM   3342  CB  LEU B 265      39.435  62.290  46.108  1.00 31.58           C  
ATOM   3343  CG  LEU B 265      39.103  61.065  45.243  1.00 31.89           C  
ATOM   3344  CD1 LEU B 265      37.656  60.580  45.378  1.00 32.01           C  
ATOM   3345  CD2 LEU B 265      39.424  61.365  43.800  1.00 31.19           C  
ATOM   3346  N   SER B 266      39.001  64.190  49.132  1.00 31.54           N  
ATOM   3347  CA  SER B 266      39.742  64.985  50.090  1.00 31.56           C  
ATOM   3348  C   SER B 266      39.663  64.338  51.471  1.00 31.67           C  
ATOM   3349  O   SER B 266      38.638  63.745  51.829  1.00 31.72           O  
ATOM   3350  CB  SER B 266      39.193  66.411  50.116  1.00 31.38           C  
ATOM   3351  OG  SER B 266      39.822  67.186  51.116  1.00 31.53           O  
ATOM   3352  N   LYS B 267      40.746  64.455  52.240  1.00 31.71           N  
ATOM   3353  CA  LYS B 267      40.807  63.920  53.601  1.00 31.81           C  
ATOM   3354  C   LYS B 267      39.790  64.607  54.517  1.00 31.89           C  
ATOM   3355  O   LYS B 267      39.512  64.132  55.620  1.00 32.17           O  
ATOM   3356  CB  LYS B 267      42.216  64.067  54.164  1.00 31.76           C  
ATOM   3357  N   ASP B 268      39.239  65.720  54.041  1.00 31.77           N  
ATOM   3358  CA  ASP B 268      38.294  66.537  54.799  1.00 31.75           C  
ATOM   3359  C   ASP B 268      36.837  66.195  54.493  1.00 31.23           C  
ATOM   3360  O   ASP B 268      35.922  66.775  55.073  1.00 31.07           O  
ATOM   3361  CB  ASP B 268      38.558  68.023  54.520  1.00 32.12           C  
ATOM   3362  CG  ASP B 268      39.920  68.482  55.032  1.00 33.59           C  
ATOM   3363  OD1 ASP B 268      40.026  69.648  55.479  1.00 34.77           O  
ATOM   3364  OD2 ASP B 268      40.882  67.673  55.006  1.00 34.99           O  
ATOM   3365  N   PHE B 269      36.627  65.250  53.581  1.00 30.68           N  
ATOM   3366  CA  PHE B 269      35.285  64.883  53.168  1.00 30.22           C  
ATOM   3367  C   PHE B 269      34.449  64.426  54.357  1.00 30.41           C  
ATOM   3368  O   PHE B 269      33.269  64.763  54.471  1.00 30.35           O  
ATOM   3369  CB  PHE B 269      35.328  63.791  52.100  1.00 29.74           C  
ATOM   3370  CG  PHE B 269      33.977  63.248  51.739  1.00 29.13           C  
ATOM   3371  CD1 PHE B 269      33.112  63.975  50.917  1.00 28.56           C  
ATOM   3372  CD2 PHE B 269      33.555  62.022  52.228  1.00 27.99           C  
ATOM   3373  CE1 PHE B 269      31.850  63.477  50.579  1.00 27.12           C  
ATOM   3374  CE2 PHE B 269      32.292  61.522  51.894  1.00 28.13           C  
ATOM   3375  CZ  PHE B 269      31.443  62.252  51.065  1.00 27.15           C  
ATOM   3376  N   SER B 270      35.081  63.662  55.240  1.00 30.65           N  
ATOM   3377  CA  SER B 270      34.400  63.038  56.369  1.00 30.95           C  
ATOM   3378  C   SER B 270      33.777  64.063  57.323  1.00 30.96           C  
ATOM   3379  O   SER B 270      32.725  63.806  57.914  1.00 31.00           O  
ATOM   3380  CB  SER B 270      35.371  62.126  57.112  1.00 30.86           C  
ATOM   3381  OG  SER B 270      36.525  62.855  57.476  1.00 31.68           O  
ATOM   3382  N   LYS B 271      34.420  65.222  57.450  1.00 30.92           N  
ATOM   3383  CA  LYS B 271      33.935  66.292  58.318  1.00 31.06           C  
ATOM   3384  C   LYS B 271      33.105  67.347  57.582  1.00 30.89           C  
ATOM   3385  O   LYS B 271      32.746  68.366  58.156  1.00 30.85           O  
ATOM   3386  CB  LYS B 271      35.086  66.942  59.096  1.00 31.18           C  
ATOM   3387  CG  LYS B 271      36.357  67.110  58.308  1.00 32.37           C  
ATOM   3388  CD  LYS B 271      37.252  68.174  58.905  1.00 34.77           C  
ATOM   3389  CE  LYS B 271      38.609  68.167  58.211  1.00 35.57           C  
ATOM   3390  NZ  LYS B 271      39.189  69.536  58.119  1.00 36.55           N  
HETATM 3391  N   MSE B 272      32.797  67.093  56.315  1.00 30.93           N  
HETATM 3392  CA  MSE B 272      31.843  67.912  55.578  1.00 30.92           C  
HETATM 3393  C   MSE B 272      30.440  67.714  56.140  1.00 30.11           C  
HETATM 3394  O   MSE B 272      30.154  66.659  56.698  1.00 30.08           O  
HETATM 3395  CB  MSE B 272      31.863  67.547  54.092  1.00 31.51           C  
HETATM 3396  CG  MSE B 272      33.124  67.967  53.354  1.00 34.29           C  
HETATM 3397 SE   MSE B 272      33.474  69.882  53.499  1.00 43.23          SE  
HETATM 3398  CE  MSE B 272      34.792  69.934  54.955  1.00 39.78           C  
ATOM   3399  N   PRO B 273      29.559  68.724  55.999  1.00 29.60           N  
ATOM   3400  CA  PRO B 273      28.153  68.589  56.415  1.00 29.21           C  
ATOM   3401  C   PRO B 273      27.434  67.405  55.764  1.00 28.69           C  
ATOM   3402  O   PRO B 273      27.743  67.038  54.633  1.00 28.43           O  
ATOM   3403  CB  PRO B 273      27.519  69.908  55.952  1.00 29.19           C  
ATOM   3404  CG  PRO B 273      28.650  70.869  55.916  1.00 29.39           C  
ATOM   3405  CD  PRO B 273      29.835  70.069  55.458  1.00 29.56           C  
ATOM   3406  N   ASP B 274      26.480  66.825  56.488  1.00 28.43           N  
ATOM   3407  CA  ASP B 274      25.762  65.629  56.037  1.00 28.27           C  
ATOM   3408  C   ASP B 274      25.222  65.762  54.613  1.00 28.15           C  
ATOM   3409  O   ASP B 274      25.519  64.920  53.760  1.00 28.09           O  
ATOM   3410  CB  ASP B 274      24.634  65.255  57.014  1.00 28.23           C  
ATOM   3411  CG  ASP B 274      25.144  64.566  58.283  1.00 28.39           C  
ATOM   3412  OD1 ASP B 274      26.361  64.610  58.574  1.00 28.32           O  
ATOM   3413  OD2 ASP B 274      24.311  63.976  59.002  1.00 28.84           O  
ATOM   3414  N   ASP B 275      24.459  66.824  54.350  1.00 27.78           N  
ATOM   3415  CA  ASP B 275      23.857  67.014  53.030  1.00 27.77           C  
ATOM   3416  C   ASP B 275      24.895  67.114  51.908  1.00 27.57           C  
ATOM   3417  O   ASP B 275      24.635  66.689  50.777  1.00 27.66           O  
ATOM   3418  CB  ASP B 275      22.951  68.234  53.018  1.00 27.86           C  
ATOM   3419  CG  ASP B 275      23.715  69.519  53.181  1.00 29.11           C  
ATOM   3420  OD1 ASP B 275      24.308  69.743  54.259  1.00 29.99           O  
ATOM   3421  OD2 ASP B 275      23.716  70.314  52.222  1.00 31.52           O  
ATOM   3422  N   VAL B 276      26.061  67.666  52.238  1.00 27.26           N  
ATOM   3423  CA  VAL B 276      27.157  67.855  51.289  1.00 27.04           C  
ATOM   3424  C   VAL B 276      27.735  66.515  50.860  1.00 27.12           C  
ATOM   3425  O   VAL B 276      27.997  66.294  49.678  1.00 27.15           O  
ATOM   3426  CB  VAL B 276      28.273  68.763  51.884  1.00 26.99           C  
ATOM   3427  CG1 VAL B 276      29.436  68.924  50.914  1.00 26.22           C  
ATOM   3428  CG2 VAL B 276      27.702  70.122  52.251  1.00 26.86           C  
ATOM   3429  N   GLN B 277      27.920  65.624  51.831  1.00 27.19           N  
ATOM   3430  CA  GLN B 277      28.447  64.289  51.582  1.00 26.83           C  
ATOM   3431  C   GLN B 277      27.446  63.488  50.772  1.00 27.14           C  
ATOM   3432  O   GLN B 277      27.792  62.928  49.722  1.00 27.58           O  
ATOM   3433  CB  GLN B 277      28.737  63.579  52.897  1.00 26.58           C  
ATOM   3434  CG  GLN B 277      29.793  64.237  53.755  1.00 26.18           C  
ATOM   3435  CD  GLN B 277      29.891  63.603  55.129  1.00 26.89           C  
ATOM   3436  OE1 GLN B 277      28.912  63.063  55.654  1.00 26.55           O  
ATOM   3437  NE2 GLN B 277      31.074  63.668  55.723  1.00 27.28           N  
ATOM   3438  N   THR B 278      26.206  63.444  51.262  1.00 27.05           N  
ATOM   3439  CA  THR B 278      25.093  62.842  50.538  1.00 27.06           C  
ATOM   3440  C   THR B 278      25.054  63.338  49.093  1.00 27.36           C  
ATOM   3441  O   THR B 278      24.972  62.542  48.160  1.00 27.29           O  
ATOM   3442  CB  THR B 278      23.752  63.157  51.224  1.00 27.05           C  
ATOM   3443  OG1 THR B 278      23.846  62.851  52.617  1.00 26.76           O  
ATOM   3444  CG2 THR B 278      22.617  62.359  50.610  1.00 26.33           C  
ATOM   3445  N   ARG B 279      25.135  64.654  48.913  1.00 27.85           N  
ATOM   3446  CA  ARG B 279      25.076  65.243  47.577  1.00 28.46           C  
ATOM   3447  C   ARG B 279      26.248  64.827  46.684  1.00 28.34           C  
ATOM   3448  O   ARG B 279      26.055  64.512  45.508  1.00 28.65           O  
ATOM   3449  CB  ARG B 279      25.002  66.760  47.655  1.00 28.74           C  
ATOM   3450  CG  ARG B 279      24.226  67.366  46.516  1.00 30.80           C  
ATOM   3451  CD  ARG B 279      22.836  67.813  46.963  1.00 32.62           C  
ATOM   3452  NE  ARG B 279      22.717  69.266  47.083  1.00 34.04           N  
ATOM   3453  CZ  ARG B 279      23.479  70.154  46.446  1.00 34.95           C  
ATOM   3454  NH1 ARG B 279      24.435  69.752  45.614  1.00 34.84           N  
ATOM   3455  NH2 ARG B 279      23.279  71.456  46.638  1.00 35.42           N  
ATOM   3456  N   LEU B 280      27.456  64.820  47.242  1.00 28.12           N  
ATOM   3457  CA  LEU B 280      28.647  64.524  46.459  1.00 27.71           C  
ATOM   3458  C   LEU B 280      28.626  63.081  45.997  1.00 27.72           C  
ATOM   3459  O   LEU B 280      28.875  62.798  44.818  1.00 27.95           O  
ATOM   3460  CB  LEU B 280      29.921  64.825  47.248  1.00 27.46           C  
ATOM   3461  CG  LEU B 280      30.300  66.300  47.362  1.00 27.20           C  
ATOM   3462  CD1 LEU B 280      31.367  66.487  48.427  1.00 27.10           C  
ATOM   3463  CD2 LEU B 280      30.749  66.879  46.022  1.00 26.96           C  
ATOM   3464  N   LEU B 281      28.312  62.176  46.921  1.00 27.33           N  
ATOM   3465  CA  LEU B 281      28.258  60.757  46.602  1.00 27.02           C  
ATOM   3466  C   LEU B 281      27.206  60.484  45.535  1.00 26.79           C  
ATOM   3467  O   LEU B 281      27.496  59.789  44.563  1.00 26.78           O  
ATOM   3468  CB  LEU B 281      28.021  59.897  47.855  1.00 27.02           C  
ATOM   3469  CG  LEU B 281      29.149  59.759  48.887  1.00 26.69           C  
ATOM   3470  CD1 LEU B 281      28.697  58.897  50.051  1.00 26.31           C  
ATOM   3471  CD2 LEU B 281      30.425  59.179  48.282  1.00 27.23           C  
ATOM   3472  N   ALA B 282      26.004  61.042  45.705  1.00 26.39           N  
ATOM   3473  CA  ALA B 282      24.939  60.908  44.699  1.00 26.37           C  
ATOM   3474  C   ALA B 282      25.365  61.437  43.324  1.00 26.22           C  
ATOM   3475  O   ALA B 282      25.248  60.741  42.318  1.00 26.25           O  
ATOM   3476  CB  ALA B 282      23.654  61.596  45.156  1.00 26.42           C  
ATOM   3477  N   SER B 283      25.873  62.663  43.298  1.00 25.95           N  
ATOM   3478  CA  SER B 283      26.292  63.297  42.062  1.00 25.86           C  
ATOM   3479  C   SER B 283      27.363  62.506  41.304  1.00 25.81           C  
ATOM   3480  O   SER B 283      27.359  62.484  40.072  1.00 25.88           O  
ATOM   3481  CB  SER B 283      26.798  64.709  42.345  1.00 25.90           C  
ATOM   3482  OG  SER B 283      27.214  65.341  41.150  1.00 25.85           O  
ATOM   3483  N   ILE B 284      28.264  61.857  42.037  1.00 25.42           N  
ATOM   3484  CA  ILE B 284      29.422  61.208  41.430  1.00 25.15           C  
ATOM   3485  C   ILE B 284      29.246  59.693  41.250  1.00 25.62           C  
ATOM   3486  O   ILE B 284      29.690  59.125  40.239  1.00 25.98           O  
ATOM   3487  CB  ILE B 284      30.740  61.532  42.213  1.00 25.10           C  
ATOM   3488  CG1 ILE B 284      31.043  63.033  42.158  1.00 24.41           C  
ATOM   3489  CG2 ILE B 284      31.928  60.754  41.654  1.00 24.21           C  
ATOM   3490  CD1 ILE B 284      32.093  63.487  43.154  1.00 24.80           C  
ATOM   3491  N   LEU B 285      28.599  59.044  42.215  1.00 25.64           N  
ATOM   3492  CA  LEU B 285      28.610  57.583  42.284  1.00 25.74           C  
ATOM   3493  C   LEU B 285      27.228  56.907  42.244  1.00 25.97           C  
ATOM   3494  O   LEU B 285      27.133  55.682  42.283  1.00 25.93           O  
ATOM   3495  CB  LEU B 285      29.440  57.111  43.501  1.00 25.67           C  
ATOM   3496  CG  LEU B 285      30.917  57.553  43.592  1.00 24.91           C  
ATOM   3497  CD1 LEU B 285      31.473  57.321  44.965  1.00 24.13           C  
ATOM   3498  CD2 LEU B 285      31.790  56.859  42.577  1.00 23.71           C  
ATOM   3499  N   GLN B 286      26.158  57.693  42.150  1.00 26.47           N  
ATOM   3500  CA  GLN B 286      24.814  57.118  42.101  1.00 26.48           C  
ATOM   3501  C   GLN B 286      24.255  57.258  40.700  1.00 26.68           C  
ATOM   3502  O   GLN B 286      24.339  58.326  40.098  1.00 26.74           O  
ATOM   3503  CB  GLN B 286      23.894  57.758  43.145  1.00 26.23           C  
ATOM   3504  CG  GLN B 286      22.756  56.842  43.657  1.00 26.03           C  
ATOM   3505  CD  GLN B 286      21.621  56.680  42.652  1.00 25.59           C  
ATOM   3506  OE1 GLN B 286      21.284  57.615  41.932  1.00 26.91           O  
ATOM   3507  NE2 GLN B 286      21.030  55.495  42.603  1.00 24.65           N  
ATOM   3508  N   ILE B 287      23.686  56.167  40.196  1.00 27.14           N  
ATOM   3509  CA  ILE B 287      23.264  56.053  38.785  1.00 27.37           C  
ATOM   3510  C   ILE B 287      22.217  57.079  38.320  1.00 27.90           C  
ATOM   3511  O   ILE B 287      22.227  57.464  37.153  1.00 28.39           O  
ATOM   3512  CB  ILE B 287      22.818  54.598  38.424  1.00 27.12           C  
ATOM   3513  CG1 ILE B 287      22.863  54.378  36.906  1.00 27.69           C  
ATOM   3514  CG2 ILE B 287      21.441  54.274  39.015  1.00 26.25           C  
ATOM   3515  CD1 ILE B 287      22.872  52.910  36.456  1.00 26.98           C  
ATOM   3516  N   ASP B 288      21.323  57.502  39.218  1.00 28.24           N  
ATOM   3517  CA  ASP B 288      20.291  58.505  38.900  1.00 28.65           C  
ATOM   3518  C   ASP B 288      20.621  59.821  39.577  1.00 28.88           C  
ATOM   3519  O   ASP B 288      19.792  60.738  39.595  1.00 28.74           O  
ATOM   3520  CB  ASP B 288      18.901  58.071  39.395  1.00 28.69           C  
ATOM   3521  CG  ASP B 288      18.529  56.663  38.971  1.00 29.84           C  
ATOM   3522  OD1 ASP B 288      18.320  56.447  37.754  1.00 29.98           O  
ATOM   3523  OD2 ASP B 288      18.430  55.779  39.866  1.00 30.81           O  
ATOM   3524  N   LYS B 289      21.821  59.904  40.154  1.00 29.18           N  
ATOM   3525  CA  LYS B 289      22.211  61.041  40.983  1.00 29.57           C  
ATOM   3526  C   LYS B 289      21.109  61.346  42.004  1.00 29.55           C  
ATOM   3527  O   LYS B 289      20.687  62.497  42.163  1.00 29.49           O  
ATOM   3528  CB  LYS B 289      22.498  62.280  40.122  1.00 29.47           C  
ATOM   3529  CG  LYS B 289      23.639  62.141  39.132  1.00 30.27           C  
ATOM   3530  CD  LYS B 289      23.818  63.453  38.376  1.00 31.79           C  
ATOM   3531  CE  LYS B 289      24.879  63.346  37.300  1.00 33.44           C  
ATOM   3532  NZ  LYS B 289      24.805  64.511  36.369  1.00 34.80           N  
ATOM   3533  N   ASP B 290      20.628  60.303  42.671  1.00 29.46           N  
ATOM   3534  CA  ASP B 290      19.517  60.445  43.601  1.00 29.25           C  
ATOM   3535  C   ASP B 290      20.032  60.408  45.024  1.00 29.14           C  
ATOM   3536  O   ASP B 290      20.460  59.362  45.496  1.00 29.42           O  
ATOM   3537  CB  ASP B 290      18.494  59.328  43.394  1.00 29.20           C  
ATOM   3538  CG  ASP B 290      17.197  59.578  44.134  1.00 29.45           C  
ATOM   3539  OD1 ASP B 290      17.103  60.564  44.890  1.00 31.20           O  
ATOM   3540  OD2 ASP B 290      16.256  58.783  43.965  1.00 30.16           O  
ATOM   3541  N   PRO B 291      19.979  61.549  45.725  1.00 29.20           N  
ATOM   3542  CA  PRO B 291      20.467  61.584  47.107  1.00 29.01           C  
ATOM   3543  C   PRO B 291      19.663  60.700  48.065  1.00 28.76           C  
ATOM   3544  O   PRO B 291      20.177  60.340  49.117  1.00 28.84           O  
ATOM   3545  CB  PRO B 291      20.316  63.060  47.501  1.00 29.05           C  
ATOM   3546  CG  PRO B 291      20.140  63.806  46.216  1.00 29.05           C  
ATOM   3547  CD  PRO B 291      19.464  62.859  45.288  1.00 29.15           C  
ATOM   3548  N   HIS B 292      18.423  60.360  47.708  1.00 28.44           N  
ATOM   3549  CA  HIS B 292      17.592  59.497  48.548  1.00 28.27           C  
ATOM   3550  C   HIS B 292      18.091  58.051  48.591  1.00 27.87           C  
ATOM   3551  O   HIS B 292      17.812  57.317  49.541  1.00 27.59           O  
ATOM   3552  CB  HIS B 292      16.140  59.520  48.084  1.00 28.55           C  
ATOM   3553  CG  HIS B 292      15.505  60.870  48.157  1.00 29.92           C  
ATOM   3554  ND1 HIS B 292      14.903  61.347  49.309  1.00 31.11           N  
ATOM   3555  CD2 HIS B 292      15.371  61.846  47.220  1.00 30.80           C  
ATOM   3556  CE1 HIS B 292      14.427  62.559  49.078  1.00 31.32           C  
ATOM   3557  NE2 HIS B 292      14.697  62.884  47.819  1.00 31.36           N  
ATOM   3558  N   LYS B 293      18.823  57.636  47.564  1.00 27.39           N  
ATOM   3559  CA  LYS B 293      19.311  56.268  47.532  1.00 27.16           C  
ATOM   3560  C   LYS B 293      20.695  56.086  48.164  1.00 26.88           C  
ATOM   3561  O   LYS B 293      21.245  54.990  48.161  1.00 26.63           O  
ATOM   3562  CB  LYS B 293      19.192  55.680  46.124  1.00 27.17           C  
ATOM   3563  CG  LYS B 293      17.834  54.994  45.917  1.00 28.01           C  
ATOM   3564  CD  LYS B 293      17.729  54.248  44.604  1.00 29.25           C  
ATOM   3565  CE  LYS B 293      16.952  55.050  43.568  1.00 29.51           C  
ATOM   3566  NZ  LYS B 293      17.132  54.496  42.196  1.00 29.89           N  
ATOM   3567  N   VAL B 294      21.223  57.164  48.745  1.00 26.85           N  
ATOM   3568  CA  VAL B 294      22.518  57.159  49.427  1.00 26.60           C  
ATOM   3569  C   VAL B 294      22.306  57.290  50.935  1.00 26.60           C  
ATOM   3570  O   VAL B 294      21.496  58.114  51.374  1.00 26.70           O  
ATOM   3571  CB  VAL B 294      23.414  58.326  48.938  1.00 26.73           C  
ATOM   3572  CG1 VAL B 294      24.798  58.271  49.588  1.00 26.77           C  
ATOM   3573  CG2 VAL B 294      23.545  58.300  47.423  1.00 26.71           C  
ATOM   3574  N   ASP B 295      23.032  56.480  51.714  1.00 26.17           N  
ATOM   3575  CA  ASP B 295      22.961  56.496  53.180  1.00 25.98           C  
ATOM   3576  C   ASP B 295      24.349  56.741  53.768  1.00 26.10           C  
ATOM   3577  O   ASP B 295      25.163  55.826  53.864  1.00 26.10           O  
ATOM   3578  CB  ASP B 295      22.374  55.177  53.711  1.00 25.91           C  
ATOM   3579  CG  ASP B 295      22.254  55.136  55.240  1.00 25.57           C  
ATOM   3580  OD1 ASP B 295      22.554  56.138  55.933  1.00 25.19           O  
ATOM   3581  OD2 ASP B 295      21.840  54.078  55.755  1.00 25.59           O  
ATOM   3582  N   ILE B 296      24.605  57.978  54.182  1.00 26.27           N  
ATOM   3583  CA  ILE B 296      25.944  58.376  54.637  1.00 26.27           C  
ATOM   3584  C   ILE B 296      26.302  57.865  56.037  1.00 26.16           C  
ATOM   3585  O   ILE B 296      27.467  57.913  56.443  1.00 25.87           O  
ATOM   3586  CB  ILE B 296      26.166  59.912  54.510  1.00 26.25           C  
ATOM   3587  CG1 ILE B 296      25.074  60.705  55.234  1.00 26.61           C  
ATOM   3588  CG2 ILE B 296      26.191  60.313  53.055  1.00 26.37           C  
ATOM   3589  CD1 ILE B 296      25.528  61.346  56.510  1.00 26.77           C  
ATOM   3590  N   LYS B 297      25.302  57.359  56.758  1.00 26.12           N  
ATOM   3591  CA  LYS B 297      25.520  56.796  58.090  1.00 26.20           C  
ATOM   3592  C   LYS B 297      26.194  55.422  58.026  1.00 26.19           C  
ATOM   3593  O   LYS B 297      26.412  54.779  59.046  1.00 26.17           O  
ATOM   3594  CB  LYS B 297      24.211  56.766  58.901  1.00 26.20           C  
ATOM   3595  CG  LYS B 297      23.708  58.157  59.352  1.00 26.78           C  
ATOM   3596  CD  LYS B 297      24.793  58.955  60.114  1.00 27.75           C  
ATOM   3597  CE  LYS B 297      24.633  60.486  59.985  1.00 27.27           C  
ATOM   3598  NZ  LYS B 297      23.721  61.101  60.999  1.00 26.58           N  
ATOM   3599  N   LYS B 298      26.535  54.986  56.820  1.00 26.21           N  
ATOM   3600  CA  LYS B 298      27.267  53.749  56.640  1.00 26.42           C  
ATOM   3601  C   LYS B 298      28.725  54.044  56.352  1.00 26.32           C  
ATOM   3602  O   LYS B 298      29.515  53.135  56.141  1.00 26.61           O  
ATOM   3603  CB  LYS B 298      26.665  52.924  55.504  1.00 26.58           C  
ATOM   3604  CG  LYS B 298      25.179  52.606  55.670  1.00 27.68           C  
ATOM   3605  CD  LYS B 298      24.917  51.591  56.769  1.00 29.29           C  
ATOM   3606  CE  LYS B 298      25.165  50.165  56.296  1.00 30.52           C  
ATOM   3607  NZ  LYS B 298      24.040  49.626  55.470  1.00 31.01           N  
ATOM   3608  N   LEU B 299      29.083  55.318  56.343  1.00 26.41           N  
ATOM   3609  CA  LEU B 299      30.461  55.711  56.095  1.00 26.64           C  
ATOM   3610  C   LEU B 299      31.362  55.422  57.287  1.00 26.98           C  
ATOM   3611  O   LEU B 299      30.954  55.538  58.447  1.00 26.83           O  
ATOM   3612  CB  LEU B 299      30.548  57.189  55.729  1.00 26.62           C  
ATOM   3613  CG  LEU B 299      29.984  57.648  54.373  1.00 26.97           C  
ATOM   3614  CD1 LEU B 299      29.974  59.171  54.316  1.00 26.44           C  
ATOM   3615  CD2 LEU B 299      30.744  57.060  53.174  1.00 25.75           C  
ATOM   3616  N   HIS B 300      32.597  55.037  56.995  1.00 27.41           N  
ATOM   3617  CA  HIS B 300      33.576  54.893  58.044  1.00 27.74           C  
ATOM   3618  C   HIS B 300      34.512  56.072  58.046  1.00 27.42           C  
ATOM   3619  O   HIS B 300      35.366  56.217  57.181  1.00 27.35           O  
ATOM   3620  CB  HIS B 300      34.271  53.539  57.978  1.00 28.26           C  
ATOM   3621  CG  HIS B 300      33.498  52.453  58.663  1.00 30.47           C  
ATOM   3622  ND1 HIS B 300      33.175  51.261  58.049  1.00 32.62           N  
ATOM   3623  CD2 HIS B 300      32.948  52.398  59.902  1.00 32.13           C  
ATOM   3624  CE1 HIS B 300      32.481  50.510  58.888  1.00 33.70           C  
ATOM   3625  NE2 HIS B 300      32.323  51.179  60.016  1.00 33.55           N  
ATOM   3626  N   LEU B 301      34.309  56.926  59.044  1.00 27.46           N  
ATOM   3627  CA  LEU B 301      34.803  58.296  59.025  1.00 27.41           C  
ATOM   3628  C   LEU B 301      36.311  58.439  59.223  1.00 27.37           C  
ATOM   3629  O   LEU B 301      36.880  59.484  58.909  1.00 27.41           O  
ATOM   3630  CB  LEU B 301      34.024  59.151  60.030  1.00 27.61           C  
ATOM   3631  CG  LEU B 301      32.519  59.361  59.788  1.00 27.97           C  
ATOM   3632  CD1 LEU B 301      31.839  59.831  61.068  1.00 28.21           C  
ATOM   3633  CD2 LEU B 301      32.240  60.337  58.640  1.00 27.82           C  
ATOM   3634  N   ASP B 302      36.957  57.393  59.728  1.00 27.33           N  
ATOM   3635  CA  ASP B 302      38.412  57.388  59.850  1.00 27.35           C  
ATOM   3636  C   ASP B 302      39.058  56.576  58.735  1.00 26.98           C  
ATOM   3637  O   ASP B 302      40.239  56.248  58.804  1.00 26.75           O  
ATOM   3638  CB  ASP B 302      38.845  56.870  61.225  1.00 27.62           C  
ATOM   3639  CG  ASP B 302      38.606  57.882  62.328  1.00 28.66           C  
ATOM   3640  OD1 ASP B 302      39.225  58.969  62.290  1.00 29.85           O  
ATOM   3641  OD2 ASP B 302      37.803  57.590  63.240  1.00 29.75           O  
ATOM   3642  N   GLY B 303      38.278  56.269  57.702  1.00 26.94           N  
ATOM   3643  CA  GLY B 303      38.758  55.487  56.562  1.00 26.99           C  
ATOM   3644  C   GLY B 303      39.867  56.175  55.790  1.00 27.04           C  
ATOM   3645  O   GLY B 303      39.799  57.379  55.546  1.00 27.03           O  
ATOM   3646  N   LYS B 304      40.898  55.420  55.421  1.00 27.16           N  
ATOM   3647  CA  LYS B 304      41.972  55.970  54.597  1.00 27.66           C  
ATOM   3648  C   LYS B 304      41.480  56.139  53.154  1.00 27.80           C  
ATOM   3649  O   LYS B 304      40.633  55.375  52.687  1.00 27.78           O  
ATOM   3650  CB  LYS B 304      43.227  55.083  54.644  1.00 27.59           C  
ATOM   3651  CG  LYS B 304      43.799  54.853  56.043  1.00 28.17           C  
ATOM   3652  CD  LYS B 304      44.898  53.775  56.072  1.00 28.04           C  
ATOM   3653  CE  LYS B 304      46.305  54.379  56.203  1.00 28.84           C  
ATOM   3654  NZ  LYS B 304      47.404  53.355  56.252  1.00 28.38           N  
ATOM   3655  N   LEU B 305      41.996  57.154  52.466  1.00 28.00           N  
ATOM   3656  CA  LEU B 305      41.726  57.352  51.046  1.00 28.12           C  
ATOM   3657  C   LEU B 305      43.030  57.208  50.294  1.00 28.32           C  
ATOM   3658  O   LEU B 305      43.974  57.957  50.529  1.00 28.49           O  
ATOM   3659  CB  LEU B 305      41.112  58.738  50.786  1.00 28.16           C  
ATOM   3660  CG  LEU B 305      40.596  59.144  49.392  1.00 28.10           C  
ATOM   3661  CD1 LEU B 305      41.692  59.731  48.512  1.00 28.73           C  
ATOM   3662  CD2 LEU B 305      39.908  57.982  48.660  1.00 28.26           C  
ATOM   3663  N   ARG B 306      43.079  56.231  49.398  1.00 28.73           N  
ATOM   3664  CA  ARG B 306      44.262  55.990  48.578  1.00 29.09           C  
ATOM   3665  C   ARG B 306      43.898  55.875  47.101  1.00 29.31           C  
ATOM   3666  O   ARG B 306      42.816  55.398  46.747  1.00 29.26           O  
ATOM   3667  CB  ARG B 306      44.976  54.726  49.039  1.00 29.09           C  
ATOM   3668  CG  ARG B 306      45.626  54.831  50.415  1.00 29.83           C  
ATOM   3669  CD  ARG B 306      45.673  53.476  51.120  1.00 30.34           C  
ATOM   3670  NE  ARG B 306      45.663  52.370  50.167  1.00 31.11           N  
ATOM   3671  CZ  ARG B 306      46.729  51.935  49.500  1.00 32.17           C  
ATOM   3672  NH1 ARG B 306      47.918  52.504  49.676  1.00 32.42           N  
ATOM   3673  NH2 ARG B 306      46.604  50.925  48.649  1.00 32.70           N  
ATOM   3674  N   PHE B 307      44.804  56.328  46.244  1.00 29.63           N  
ATOM   3675  CA  PHE B 307      44.598  56.230  44.808  1.00 30.09           C  
ATOM   3676  C   PHE B 307      45.905  56.077  44.043  1.00 30.24           C  
ATOM   3677  O   PHE B 307      46.972  56.461  44.528  1.00 30.09           O  
ATOM   3678  CB  PHE B 307      43.779  57.413  44.279  1.00 30.12           C  
ATOM   3679  CG  PHE B 307      44.417  58.747  44.500  1.00 30.20           C  
ATOM   3680  CD1 PHE B 307      45.316  59.264  43.575  1.00 30.57           C  
ATOM   3681  CD2 PHE B 307      44.101  59.503  45.618  1.00 30.71           C  
ATOM   3682  CE1 PHE B 307      45.905  60.511  43.771  1.00 30.37           C  
ATOM   3683  CE2 PHE B 307      44.687  60.755  45.821  1.00 30.83           C  
ATOM   3684  CZ  PHE B 307      45.591  61.255  44.895  1.00 30.24           C  
ATOM   3685  N   ALA B 308      45.806  55.507  42.848  1.00 30.52           N  
ATOM   3686  CA  ALA B 308      46.967  55.320  41.987  1.00 31.11           C  
ATOM   3687  C   ALA B 308      46.626  55.569  40.521  1.00 31.29           C  
ATOM   3688  O   ALA B 308      45.601  55.106  40.016  1.00 31.58           O  
ATOM   3689  CB  ALA B 308      47.541  53.924  42.170  1.00 31.29           C  
ATOM   3690  N   SER B 309      47.493  56.300  39.840  1.00 31.51           N  
ATOM   3691  CA  SER B 309      47.271  56.616  38.443  1.00 31.86           C  
ATOM   3692  C   SER B 309      48.210  55.835  37.537  1.00 32.18           C  
ATOM   3693  O   SER B 309      49.232  55.302  37.978  1.00 32.27           O  
ATOM   3694  CB  SER B 309      47.451  58.107  38.204  1.00 31.52           C  
ATOM   3695  OG  SER B 309      48.795  58.464  38.429  1.00 32.00           O  
ATOM   3696  N   HIS B 310      47.838  55.778  36.264  1.00 32.58           N  
ATOM   3697  CA  HIS B 310      48.637  55.162  35.228  1.00 33.03           C  
ATOM   3698  C   HIS B 310      48.152  55.700  33.888  1.00 33.48           C  
ATOM   3699  O   HIS B 310      46.937  55.802  33.660  1.00 33.64           O  
ATOM   3700  CB  HIS B 310      48.478  53.641  35.265  1.00 32.80           C  
ATOM   3701  CG  HIS B 310      49.600  52.901  34.608  1.00 33.18           C  
ATOM   3702  ND1 HIS B 310      49.833  52.955  33.249  1.00 33.55           N  
ATOM   3703  CD2 HIS B 310      50.556  52.091  35.123  1.00 33.64           C  
ATOM   3704  CE1 HIS B 310      50.889  52.217  32.956  1.00 34.01           C  
ATOM   3705  NE2 HIS B 310      51.345  51.680  34.075  1.00 34.86           N  
ATOM   3706  N   GLU B 311      49.094  56.053  33.013  1.00 33.64           N  
ATOM   3707  CA  GLU B 311      48.772  56.451  31.639  1.00 33.93           C  
ATOM   3708  C   GLU B 311      48.387  55.251  30.773  1.00 33.72           C  
ATOM   3709  O   GLU B 311      48.904  54.151  30.963  1.00 33.58           O  
ATOM   3710  CB  GLU B 311      49.945  57.178  30.999  1.00 34.05           C  
ATOM   3711  CG  GLU B 311      51.219  56.352  30.954  1.00 36.01           C  
ATOM   3712  CD  GLU B 311      52.325  57.001  30.147  1.00 38.73           C  
ATOM   3713  OE1 GLU B 311      52.163  58.163  29.709  1.00 39.39           O  
ATOM   3714  OE2 GLU B 311      53.366  56.336  29.949  1.00 40.33           O  
ATOM   3715  N   TYR B 312      47.476  55.471  29.829  1.00 33.53           N  
ATOM   3716  CA  TYR B 312      47.114  54.435  28.871  1.00 33.47           C  
ATOM   3717  C   TYR B 312      47.034  55.008  27.463  1.00 33.50           C  
ATOM   3718  O   TYR B 312      46.729  56.189  27.286  1.00 33.39           O  
ATOM   3719  CB  TYR B 312      45.803  53.727  29.274  1.00 33.41           C  
ATOM   3720  CG  TYR B 312      44.543  54.566  29.182  1.00 33.21           C  
ATOM   3721  CD1 TYR B 312      43.740  54.537  28.036  1.00 32.89           C  
ATOM   3722  CD2 TYR B 312      44.139  55.370  30.247  1.00 33.23           C  
ATOM   3723  CE1 TYR B 312      42.579  55.301  27.947  1.00 32.33           C  
ATOM   3724  CE2 TYR B 312      42.977  56.141  30.172  1.00 32.82           C  
ATOM   3725  CZ  TYR B 312      42.206  56.105  29.018  1.00 33.13           C  
ATOM   3726  OH  TYR B 312      41.059  56.865  28.942  1.00 33.19           O  
ATOM   3727  N   ASP B 313      47.345  54.182  26.466  1.00 33.51           N  
ATOM   3728  CA  ASP B 313      47.086  54.552  25.083  1.00 33.58           C  
ATOM   3729  C   ASP B 313      45.688  54.087  24.717  1.00 33.84           C  
ATOM   3730  O   ASP B 313      45.302  52.955  25.034  1.00 33.83           O  
ATOM   3731  CB  ASP B 313      48.097  53.932  24.126  1.00 33.61           C  
ATOM   3732  CG  ASP B 313      47.782  54.256  22.671  1.00 34.18           C  
ATOM   3733  OD1 ASP B 313      48.082  55.391  22.220  1.00 34.17           O  
ATOM   3734  OD2 ASP B 313      47.212  53.381  21.988  1.00 33.78           O  
ATOM   3735  N   PHE B 314      44.928  54.956  24.050  1.00 33.90           N  
ATOM   3736  CA  PHE B 314      43.552  54.627  23.684  1.00 33.91           C  
ATOM   3737  C   PHE B 314      43.447  53.387  22.779  1.00 33.55           C  
ATOM   3738  O   PHE B 314      42.616  52.511  23.019  1.00 33.50           O  
ATOM   3739  CB  PHE B 314      42.837  55.838  23.078  1.00 34.19           C  
ATOM   3740  CG  PHE B 314      41.350  55.658  22.940  1.00 34.89           C  
ATOM   3741  CD1 PHE B 314      40.569  55.307  24.041  1.00 35.38           C  
ATOM   3742  CD2 PHE B 314      40.724  55.853  21.711  1.00 35.57           C  
ATOM   3743  CE1 PHE B 314      39.185  55.144  23.918  1.00 35.98           C  
ATOM   3744  CE2 PHE B 314      39.344  55.694  21.578  1.00 36.23           C  
ATOM   3745  CZ  PHE B 314      38.574  55.335  22.680  1.00 35.60           C  
ATOM   3746  N   ARG B 315      44.301  53.300  21.764  1.00 33.33           N  
ATOM   3747  CA  ARG B 315      44.289  52.146  20.862  1.00 33.31           C  
ATOM   3748  C   ARG B 315      44.559  50.835  21.577  1.00 32.83           C  
ATOM   3749  O   ARG B 315      43.845  49.857  21.366  1.00 32.81           O  
ATOM   3750  CB  ARG B 315      45.286  52.316  19.719  1.00 33.24           C  
ATOM   3751  CG  ARG B 315      44.677  52.922  18.510  1.00 34.47           C  
ATOM   3752  CD  ARG B 315      43.879  51.912  17.678  1.00 36.28           C  
ATOM   3753  NE  ARG B 315      43.474  52.545  16.429  1.00 36.12           N  
ATOM   3754  CZ  ARG B 315      44.282  52.708  15.392  1.00 37.78           C  
ATOM   3755  NH1 ARG B 315      45.530  52.252  15.446  1.00 38.24           N  
ATOM   3756  NH2 ARG B 315      43.845  53.321  14.295  1.00 39.29           N  
ATOM   3757  N   GLN B 316      45.589  50.815  22.414  1.00 32.28           N  
ATOM   3758  CA  GLN B 316      45.928  49.606  23.143  1.00 32.08           C  
ATOM   3759  C   GLN B 316      44.718  49.103  23.945  1.00 31.99           C  
ATOM   3760  O   GLN B 316      44.366  47.920  23.875  1.00 32.14           O  
ATOM   3761  CB  GLN B 316      47.147  49.836  24.032  1.00 31.91           C  
ATOM   3762  CG  GLN B 316      47.744  48.562  24.590  1.00 31.99           C  
ATOM   3763  CD  GLN B 316      47.990  47.521  23.530  1.00 32.31           C  
ATOM   3764  OE1 GLN B 316      47.304  46.502  23.481  1.00 33.38           O  
ATOM   3765  NE2 GLN B 316      48.961  47.772  22.662  1.00 32.11           N  
ATOM   3766  N   PHE B 317      44.071  50.017  24.668  1.00 31.79           N  
ATOM   3767  CA  PHE B 317      42.837  49.731  25.389  1.00 31.67           C  
ATOM   3768  C   PHE B 317      41.781  49.134  24.460  1.00 31.83           C  
ATOM   3769  O   PHE B 317      41.169  48.116  24.796  1.00 31.68           O  
ATOM   3770  CB  PHE B 317      42.315  50.994  26.082  1.00 31.48           C  
ATOM   3771  CG  PHE B 317      40.860  50.925  26.478  1.00 31.49           C  
ATOM   3772  CD1 PHE B 317      40.394  49.933  27.325  1.00 31.38           C  
ATOM   3773  CD2 PHE B 317      39.960  51.866  26.009  1.00 32.26           C  
ATOM   3774  CE1 PHE B 317      39.054  49.871  27.681  1.00 31.57           C  
ATOM   3775  CE2 PHE B 317      38.618  51.807  26.372  1.00 32.50           C  
ATOM   3776  CZ  PHE B 317      38.174  50.806  27.215  1.00 31.08           C  
ATOM   3777  N   GLN B 318      41.582  49.762  23.298  1.00 31.91           N  
ATOM   3778  CA  GLN B 318      40.682  49.239  22.265  1.00 32.02           C  
ATOM   3779  C   GLN B 318      41.043  47.818  21.835  1.00 32.00           C  
ATOM   3780  O   GLN B 318      40.156  47.020  21.535  1.00 32.26           O  
ATOM   3781  CB  GLN B 318      40.718  50.118  21.029  1.00 32.08           C  
ATOM   3782  CG  GLN B 318      39.842  51.348  21.062  1.00 33.10           C  
ATOM   3783  CD  GLN B 318      40.053  52.184  19.816  1.00 33.86           C  
ATOM   3784  OE1 GLN B 318      41.079  52.849  19.670  1.00 34.29           O  
ATOM   3785  NE2 GLN B 318      39.098  52.127  18.895  1.00 33.82           N  
ATOM   3786  N   ARG B 319      42.338  47.514  21.788  1.00 31.72           N  
ATOM   3787  CA  ARG B 319      42.797  46.206  21.347  1.00 31.68           C  
ATOM   3788  C   ARG B 319      42.461  45.176  22.402  1.00 31.63           C  
ATOM   3789  O   ARG B 319      41.924  44.105  22.101  1.00 31.65           O  
ATOM   3790  CB  ARG B 319      44.306  46.209  21.097  1.00 31.79           C  
ATOM   3791  CG  ARG B 319      44.748  46.900  19.816  1.00 32.15           C  
ATOM   3792  CD  ARG B 319      46.217  46.605  19.562  1.00 33.04           C  
ATOM   3793  NE  ARG B 319      46.781  47.537  18.597  1.00 33.75           N  
ATOM   3794  CZ  ARG B 319      47.653  48.500  18.886  1.00 34.00           C  
ATOM   3795  NH1 ARG B 319      48.098  48.663  20.124  1.00 33.33           N  
ATOM   3796  NH2 ARG B 319      48.085  49.303  17.919  1.00 33.79           N  
ATOM   3797  N   ASN B 320      42.788  45.518  23.644  1.00 31.33           N  
ATOM   3798  CA  ASN B 320      42.554  44.650  24.769  1.00 30.75           C  
ATOM   3799  C   ASN B 320      41.057  44.389  24.947  1.00 30.82           C  
ATOM   3800  O   ASN B 320      40.642  43.246  25.130  1.00 30.81           O  
ATOM   3801  CB  ASN B 320      43.149  45.269  26.035  1.00 30.87           C  
ATOM   3802  CG  ASN B 320      44.669  45.400  25.988  1.00 29.96           C  
ATOM   3803  OD1 ASN B 320      45.328  44.916  25.067  1.00 29.71           O  
ATOM   3804  ND2 ASN B 320      45.231  46.064  26.999  1.00 28.02           N  
ATOM   3805  N   ALA B 321      40.249  45.444  24.878  1.00 30.72           N  
ATOM   3806  CA  ALA B 321      38.801  45.312  24.985  1.00 30.98           C  
ATOM   3807  C   ALA B 321      38.254  44.404  23.883  1.00 31.62           C  
ATOM   3808  O   ALA B 321      37.369  43.578  24.126  1.00 31.31           O  
ATOM   3809  CB  ALA B 321      38.149  46.666  24.917  1.00 30.94           C  
ATOM   3810  N   GLN B 322      38.811  44.564  22.682  1.00 32.21           N  
ATOM   3811  CA  GLN B 322      38.413  43.817  21.500  1.00 32.91           C  
ATOM   3812  C   GLN B 322      38.713  42.329  21.622  1.00 32.93           C  
ATOM   3813  O   GLN B 322      37.854  41.493  21.303  1.00 33.01           O  
ATOM   3814  CB  GLN B 322      39.126  44.376  20.274  1.00 33.32           C  
ATOM   3815  CG  GLN B 322      38.691  43.758  18.979  1.00 35.92           C  
ATOM   3816  CD  GLN B 322      37.339  44.253  18.542  1.00 40.64           C  
ATOM   3817  OE1 GLN B 322      36.304  43.855  19.100  1.00 42.88           O  
ATOM   3818  NE2 GLN B 322      37.328  45.139  17.534  1.00 41.27           N  
ATOM   3819  N   TYR B 323      39.927  41.998  22.070  1.00 32.76           N  
ATOM   3820  CA  TYR B 323      40.310  40.595  22.259  1.00 32.60           C  
ATOM   3821  C   TYR B 323      39.313  39.877  23.171  1.00 32.04           C  
ATOM   3822  O   TYR B 323      38.817  38.804  22.828  1.00 31.98           O  
ATOM   3823  CB  TYR B 323      41.749  40.451  22.786  1.00 32.95           C  
ATOM   3824  CG  TYR B 323      42.040  39.069  23.352  1.00 33.64           C  
ATOM   3825  CD1 TYR B 323      42.438  38.016  22.523  1.00 33.09           C  
ATOM   3826  CD2 TYR B 323      41.883  38.811  24.721  1.00 34.94           C  
ATOM   3827  CE1 TYR B 323      42.678  36.747  23.040  1.00 33.39           C  
ATOM   3828  CE2 TYR B 323      42.123  37.539  25.253  1.00 34.87           C  
ATOM   3829  CZ  TYR B 323      42.518  36.513  24.409  1.00 34.29           C  
ATOM   3830  OH  TYR B 323      42.752  35.261  24.951  1.00 34.01           O  
ATOM   3831  N   VAL B 324      39.022  40.489  24.316  1.00 31.60           N  
ATOM   3832  CA  VAL B 324      38.084  39.933  25.287  1.00 31.35           C  
ATOM   3833  C   VAL B 324      36.682  39.820  24.699  1.00 31.33           C  
ATOM   3834  O   VAL B 324      36.023  38.788  24.849  1.00 31.36           O  
ATOM   3835  CB  VAL B 324      38.014  40.765  26.589  1.00 31.05           C  
ATOM   3836  CG1 VAL B 324      36.996  40.159  27.540  1.00 30.68           C  
ATOM   3837  CG2 VAL B 324      39.378  40.849  27.246  1.00 30.18           C  
ATOM   3838  N   ALA B 325      36.244  40.877  24.023  1.00 31.04           N  
ATOM   3839  CA  ALA B 325      34.897  40.921  23.459  1.00 31.17           C  
ATOM   3840  C   ALA B 325      34.718  39.984  22.264  1.00 31.33           C  
ATOM   3841  O   ALA B 325      33.590  39.701  21.852  1.00 31.44           O  
ATOM   3842  CB  ALA B 325      34.517  42.349  23.084  1.00 30.99           C  
ATOM   3843  N   GLY B 326      35.829  39.509  21.710  1.00 31.40           N  
ATOM   3844  CA  GLY B 326      35.792  38.619  20.556  1.00 31.52           C  
ATOM   3845  C   GLY B 326      36.054  37.166  20.903  1.00 31.59           C  
ATOM   3846  O   GLY B 326      36.543  36.410  20.074  1.00 31.79           O  
ATOM   3847  N   LEU B 327      35.731  36.768  22.128  1.00 31.74           N  
ATOM   3848  CA  LEU B 327      35.869  35.375  22.534  1.00 31.36           C  
ATOM   3849  C   LEU B 327      34.498  34.710  22.589  1.00 31.55           C  
ATOM   3850  O   LEU B 327      34.233  33.754  21.860  1.00 31.68           O  
ATOM   3851  CB  LEU B 327      36.556  35.290  23.891  1.00 31.16           C  
ATOM   3852  CG  LEU B 327      38.021  35.718  23.953  1.00 30.51           C  
ATOM   3853  CD1 LEU B 327      38.328  36.260  25.333  1.00 30.14           C  
ATOM   3854  CD2 LEU B 327      38.957  34.577  23.610  1.00 29.07           C  
TER    3855      LEU B 327                                                      
HETATM 3856  O   HOH A 329      44.183  24.617  43.336  1.00 34.44           O  
HETATM 3857  O   HOH A 330      45.062  30.510  52.191  1.00 27.52           O  
HETATM 3858  O   HOH A 331      54.384  37.265  71.833  1.00 27.22           O  
HETATM 3859  O   HOH A 332      48.376  19.251  62.910  1.00 29.82           O  
HETATM 3860  O   HOH A 333      33.303  25.871  59.650  1.00 40.66           O  
HETATM 3861  O   HOH A 334      33.448  30.788  48.718  1.00 29.84           O  
HETATM 3862  O   HOH A 335      33.123  37.640  71.965  1.00 46.12           O  
HETATM 3863  O   HOH A 336      49.302  50.420  70.874  1.00 49.77           O  
HETATM 3864  O   HOH A 337      32.084  30.399  50.762  1.00 29.57           O  
HETATM 3865  O   HOH A 338      40.908  29.235  46.547  1.00 24.96           O  
HETATM 3866  O   HOH A 339      13.528  38.377  44.481  1.00 38.94           O  
HETATM 3867  O   HOH A 340      67.133  26.776  65.249  1.00 48.50           O  
HETATM 3868  O   HOH A 341      50.344  24.550  67.557  1.00 43.34           O  
HETATM 3869  O   HOH A 342      28.730  32.962  65.899  1.00 35.75           O  
HETATM 3870  O   HOH A 343      43.775  29.089  44.244  1.00 45.61           O  
HETATM 3871  O   HOH A 344      44.223  34.964  45.610  1.00 45.59           O  
HETATM 3872  O   HOH A 345      51.926  23.425  59.013  1.00 32.65           O  
HETATM 3873  O   HOH A 346      41.492  34.183  46.804  1.00 31.58           O  
HETATM 3874  O   HOH A 347      45.572  36.889  75.670  1.00 30.88           O  
HETATM 3875  O   HOH A 348      21.249  42.824  62.183  1.00 41.89           O  
HETATM 3876  O   HOH A 349      58.788  19.510  56.374  1.00 31.52           O  
HETATM 3877  O   HOH A 350      41.366  55.065  18.021  1.00 33.24           O  
HETATM 3878  O   HOH A 351      50.582  56.179  40.818  1.00 50.16           O  
HETATM 3879  O   HOH A 352      47.329  47.102  65.233  1.00 38.15           O  
HETATM 3880  O   HOH A 353      33.207  55.505  61.196  1.00 56.19           O  
HETATM 3881  O   HOH A 354      38.855  19.929  53.232  1.00 36.44           O  
HETATM 3882  O   HOH A 355      35.561  26.664  55.333  1.00 39.16           O  
HETATM 3883  O   HOH A 356      56.519  39.088  68.645  1.00 46.63           O  
HETATM 3884  O   HOH A 357      34.336  29.516  27.989  1.00 50.27           O  
HETATM 3885  O   HOH A 358      53.787  22.013  72.967  1.00 43.04           O  
HETATM 3886  O   HOH A 359      31.222  38.641  56.284  1.00 39.75           O  
HETATM 3887  O   HOH A 360      21.702  63.701  54.686  1.00 46.46           O  
HETATM 3888  O   HOH A 361      51.399  58.165  38.233  1.00 58.18           O  
HETATM 3889  O   HOH A 362      40.613  72.849  56.681  1.00 30.77           O  
HETATM 3890  O   HOH A 363      29.501  41.044  33.501  1.00 31.62           O  
HETATM 3891  O   HOH A 364      59.751  17.402  49.295  1.00 43.79           O  
HETATM 3892  O   HOH A 365      32.754  45.244  58.104  1.00 42.86           O  
HETATM 3893  O   HOH A 366      25.780  53.237  43.469  1.00 30.62           O  
HETATM 3894  O   HOH A 367      51.626  43.112  65.394  1.00 40.47           O  
HETATM 3895  O   HOH A 368      43.520  33.047  34.330  1.00 39.02           O  
HETATM 3896  O   HOH A 369      37.188  64.731  18.895  1.00 49.80           O  
HETATM 3897  O   HOH A 370      38.517  64.375  25.488  1.00 41.56           O  
HETATM 3898  O   HOH A 371      47.361  48.791  56.442  1.00 47.91           O  
HETATM 3899  O   HOH A 372      43.926  38.429  43.931  1.00 50.27           O  
HETATM 3900  O   HOH A 373      44.660  31.193  31.511  1.00 43.14           O  
HETATM 3901  O   HOH A 374      40.438  17.812  52.627  1.00 30.25           O  
HETATM 3902  O   HOH A 375      28.724  59.480  58.823  1.00 46.07           O  
HETATM 3903  O   HOH A 376      28.844  37.031  73.003  1.00 54.73           O  
HETATM 3904  O   HOH A 377      36.691  19.532  51.480  1.00 47.06           O  
HETATM 3905  O   HOH A 378      57.945  27.980  69.547  1.00 49.44           O  
HETATM 3906  O   HOH A 379      49.924  41.369  67.587  1.00 51.34           O  
HETATM 3907  O   HOH B   1      30.773  45.291  45.825  1.00 15.37           O  
HETATM 3908  O   HOH B   3      28.305  64.197  38.434  1.00 45.13           O  
HETATM 3909  O   HOH B   4      25.435  32.217  34.237  1.00 36.90           O  
HETATM 3910  O   HOH B   5      38.752  40.913  18.391  1.00 40.31           O  
HETATM 3911  O   HOH B   7      48.199  51.323  27.416  1.00 29.39           O  
HETATM 3912  O   HOH B   8      30.616  49.783  48.725  1.00 28.20           O  
HETATM 3913  O   HOH B   9      45.925  57.270  22.823  1.00 31.76           O  
HETATM 3914  O   HOH B  10      27.408  50.044  48.893  1.00 33.37           O  
HETATM 3915  O   HOH B  12      34.418  49.437  41.276  1.00 25.66           O  
HETATM 3916  O   HOH B  13      32.156  34.939  25.017  1.00 39.76           O  
HETATM 3917  O   HOH B  16      28.582  60.429  21.832  1.00 35.31           O  
HETATM 3918  O   HOH B  18      52.173  45.603  38.591  1.00 46.71           O  
HETATM 3919  O   HOH B  19      34.180  13.888  41.939  1.00 37.96           O  
HETATM 3920  O   HOH B  21      31.025  59.118  34.692  1.00 30.23           O  
HETATM 3921  O   HOH B  22      48.157  58.528  26.822  1.00 35.21           O  
HETATM 3922  O   HOH B  25      23.333  49.489  49.578  1.00 33.89           O  
HETATM 3923  O   HOH B  26      36.706  46.127  47.965  1.00 36.96           O  
HETATM 3924  O   HOH B  27      30.917  47.482  49.965  1.00 53.22           O  
HETATM 3925  O   HOH B  29      36.841  53.608  59.379  1.00 43.88           O  
CONECT  106  112                                                                
CONECT  112  106  113                                                           
CONECT  113  112  114  116                                                      
CONECT  114  113  115  120                                                      
CONECT  115  114                                                                
CONECT  116  113  117                                                           
CONECT  117  116  118                                                           
CONECT  118  117  119                                                           
CONECT  119  118                                                                
CONECT  120  114                                                                
CONECT  195  200                                                                
CONECT  200  195  201                                                           
CONECT  201  200  202  204                                                      
CONECT  202  201  203  208                                                      
CONECT  203  202                                                                
CONECT  204  201  205                                                           
CONECT  205  204  206                                                           
CONECT  206  205  207                                                           
CONECT  207  206                                                                
CONECT  208  202                                                                
CONECT  570  576                                                                
CONECT  576  570  577                                                           
CONECT  577  576  578  580                                                      
CONECT  578  577  579  584                                                      
CONECT  579  578                                                                
CONECT  580  577  581                                                           
CONECT  581  580  582                                                           
CONECT  582  581  583                                                           
CONECT  583  582                                                                
CONECT  584  578                                                                
CONECT  600  605                                                                
CONECT  605  600  606                                                           
CONECT  606  605  607  609                                                      
CONECT  607  606  608                                                           
CONECT  608  607  611  612  613                                                 
CONECT  609  606  610  614                                                      
CONECT  610  609                                                                
CONECT  611  608                                                                
CONECT  612  608                                                                
CONECT  613  608                                                                
CONECT  614  609                                                                
CONECT  683  689                                                                
CONECT  689  683  690                                                           
CONECT  690  689  691  693                                                      
CONECT  691  690  692  697                                                      
CONECT  692  691                                                                
CONECT  693  690  694                                                           
CONECT  694  693  695                                                           
CONECT  695  694  696                                                           
CONECT  696  695                                                                
CONECT  697  691                                                                
CONECT 1468 1475                                                                
CONECT 1475 1468 1476                                                           
CONECT 1476 1475 1477 1479                                                      
CONECT 1477 1476 1478 1483                                                      
CONECT 1478 1477                                                                
CONECT 1479 1476 1480                                                           
CONECT 1480 1479 1481                                                           
CONECT 1481 1480 1482                                                           
CONECT 1482 1481                                                                
CONECT 1483 1477                                                                
CONECT 2025 2031                                                                
CONECT 2031 2025 2032                                                           
CONECT 2032 2031 2033 2035                                                      
CONECT 2033 2032 2034 2039                                                      
CONECT 2034 2033                                                                
CONECT 2035 2032 2036                                                           
CONECT 2036 2035 2037                                                           
CONECT 2037 2036 2038                                                           
CONECT 2038 2037                                                                
CONECT 2039 2033                                                                
CONECT 2114 2119                                                                
CONECT 2119 2114 2120                                                           
CONECT 2120 2119 2121 2123                                                      
CONECT 2121 2120 2122 2127                                                      
CONECT 2122 2121                                                                
CONECT 2123 2120 2124                                                           
CONECT 2124 2123 2125                                                           
CONECT 2125 2124 2126                                                           
CONECT 2126 2125                                                                
CONECT 2127 2121                                                                
CONECT 2483 2489                                                                
CONECT 2489 2483 2490                                                           
CONECT 2490 2489 2491 2493                                                      
CONECT 2491 2490 2492 2497                                                      
CONECT 2492 2491                                                                
CONECT 2493 2490 2494                                                           
CONECT 2494 2493 2495                                                           
CONECT 2495 2494 2496                                                           
CONECT 2496 2495                                                                
CONECT 2497 2491                                                                
CONECT 2513 2518                                                                
CONECT 2518 2513 2519                                                           
CONECT 2519 2518 2520 2522                                                      
CONECT 2520 2519 2521                                                           
CONECT 2521 2520 2524 2525 2526                                                 
CONECT 2522 2519 2523 2527                                                      
CONECT 2523 2522                                                                
CONECT 2524 2521                                                                
CONECT 2525 2521                                                                
CONECT 2526 2521                                                                
CONECT 2527 2522                                                                
CONECT 2596 2602                                                                
CONECT 2602 2596 2603                                                           
CONECT 2603 2602 2604 2606                                                      
CONECT 2604 2603 2605 2610                                                      
CONECT 2605 2604                                                                
CONECT 2606 2603 2607                                                           
CONECT 2607 2606 2608                                                           
CONECT 2608 2607 2609                                                           
CONECT 2609 2608                                                                
CONECT 2610 2604                                                                
CONECT 3384 3391                                                                
CONECT 3391 3384 3392                                                           
CONECT 3392 3391 3393 3395                                                      
CONECT 3393 3392 3394 3399                                                      
CONECT 3394 3393                                                                
CONECT 3395 3392 3396                                                           
CONECT 3396 3395 3397                                                           
CONECT 3397 3396 3398                                                           
CONECT 3398 3397                                                                
CONECT 3399 3393                                                                
MASTER      605    0   12   29    8    0    0    6 3923    2  122   44          
END                                                                             

    </textarea>
  
*/
