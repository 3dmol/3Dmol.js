
  /* @div
  <div style="width: 400px; height: 400px; position: relative;"
               class='viewer_3Dmoljs'
               data-element='moldata'
               data-backgroundcolor='0xffffff'
               data-style='cartoon:color=spectrum,colorscheme=roygb'></div>

               */

               /*
 
@data <textarea style="display: none;" id="moldata">
COMPND    CRAMBIN
ATOM      1  N   THR A   1      17.047  14.099   3.625  1.00 13.79           N  
ATOM      2  CA  THR A   1      16.967  12.784   4.338  1.00 10.80           C  
ATOM      3  C   THR A   1      15.685  12.755   5.133  1.00  9.19           C  
ATOM      4  O   THR A   1      15.268  13.825   5.594  1.00  9.85           O  
ATOM      5  CB  THR A   1      18.170  12.703   5.337  1.00 13.02           C  
ATOM      6  OG1 THR A   1      19.334  12.829   4.463  1.00 15.06           O  
ATOM      7  CG2 THR A   1      18.150  11.546   6.304  1.00 14.23           C  
ATOM      8  N   THR A   2      15.115  11.555   5.265  1.00  7.81           N  
ATOM      9  CA  THR A   2      13.856  11.469   6.066  1.00  8.31           C  
ATOM     10  C   THR A   2      14.164  10.785   7.379  1.00  5.80           C  
ATOM     11  O   THR A   2      14.993   9.862   7.443  1.00  6.94           O  
ATOM     12  CB  THR A   2      12.732  10.711   5.261  1.00 10.32           C  
ATOM     13  OG1 THR A   2      13.308   9.439   4.926  1.00 12.81           O  
ATOM     14  CG2 THR A   2      12.484  11.442   3.895  1.00 11.90           C  
ATOM     15  N   CYS A   3      13.488  11.241   8.417  1.00  5.24           N  
ATOM     16  CA  CYS A   3      13.660  10.707   9.787  1.00  5.39           C  
ATOM     17  C   CYS A   3      12.269  10.431  10.323  1.00  4.45           C  
ATOM     18  O   CYS A   3      11.393  11.308  10.185  1.00  6.54           O  
ATOM     19  CB  CYS A   3      14.368  11.748  10.691  1.00  5.99           C  
ATOM     20  SG  CYS A   3      15.885  12.426  10.016  1.00  7.01           S  
ATOM     21  N   CYS A   4      12.019   9.272  10.928  1.00  3.90           N  
ATOM     22  CA  CYS A   4      10.646   8.991  11.408  1.00  4.24           C  
ATOM     23  C   CYS A   4      10.654   8.793  12.919  1.00  3.72           C  
ATOM     24  O   CYS A   4      11.659   8.296  13.491  1.00  5.30           O  
ATOM     25  CB  CYS A   4      10.057   7.752  10.682  1.00  4.41           C  
ATOM     26  SG  CYS A   4       9.837   8.018   8.904  1.00  4.72           S  
ATOM     27  N   PRO A   5       9.561   9.108  13.563  1.00  3.96           N  
ATOM     28  CA  PRO A   5       9.448   9.034  15.012  1.00  4.25           C  
ATOM     29  C   PRO A   5       9.288   7.670  15.606  1.00  4.96           C  
ATOM     30  O   PRO A   5       9.490   7.519  16.819  1.00  7.44           O  
ATOM     31  CB  PRO A   5       8.230   9.957  15.345  1.00  5.11           C  
ATOM     32  CG  PRO A   5       7.338   9.786  14.114  1.00  5.24           C  
ATOM     33  CD  PRO A   5       8.366   9.804  12.958  1.00  5.20           C  
ATOM     34  N   SER A   6       8.875   6.686  14.796  1.00  4.83           N  
ATOM     35  CA  SER A   6       8.673   5.314  15.279  1.00  4.45           C  
ATOM     36  C   SER A   6       8.753   4.376  14.083  1.00  4.99           C  
ATOM     37  O   SER A   6       8.726   4.858  12.923  1.00  4.61           O  
ATOM     38  CB  SER A   6       7.340   5.121  15.996  1.00  5.05           C  
ATOM     39  OG  SER A   6       6.274   5.220  15.031  1.00  6.39           O  
ATOM     40  N   ILE A   7       8.881   3.075  14.358  1.00  4.94           N  
ATOM     41  CA  ILE A   7       8.912   2.083  13.258  1.00  6.33           C  
ATOM     42  C   ILE A   7       7.581   2.090  12.506  1.00  5.32           C  
ATOM     43  O   ILE A   7       7.670   2.031  11.245  1.00  6.85           O  
ATOM     44  CB  ILE A   7       9.207   0.677  13.924  1.00  8.43           C  
ATOM     45  CG1 ILE A   7      10.714   0.702  14.312  1.00  9.78           C  
ATOM     46  CG2 ILE A   7       8.811  -0.477  12.969  1.00 11.70           C  
ATOM     47  CD1 ILE A   7      11.185  -0.516  15.142  1.00  9.92           C  
ATOM     48  N   VAL A   8       6.458   2.162  13.159  1.00  5.02           N  
ATOM     49  CA  VAL A   8       5.145   2.209  12.453  1.00  6.93           C  
ATOM     50  C   VAL A   8       5.115   3.379  11.461  1.00  5.39           C  
ATOM     51  O   VAL A   8       4.664   3.268  10.343  1.00  6.30           O  
ATOM     52  CB  VAL A   8       3.995   2.354  13.478  1.00  9.64           C  
ATOM     53  CG1 VAL A   8       2.716   2.891  12.869  1.00 13.85           C  
ATOM     54  CG2 VAL A   8       3.758   1.032  14.208  1.00 11.97           C  
ATOM     55  N   ALA A   9       5.606   4.546  11.941  1.00  3.73           N  
ATOM     56  CA  ALA A   9       5.598   5.767  11.082  1.00  3.56           C  
ATOM     57  C   ALA A   9       6.441   5.527   9.850  1.00  4.13           C  
ATOM     58  O   ALA A   9       6.052   5.933   8.744  1.00  4.36           O  
ATOM     59  CB  ALA A   9       6.022   6.977  11.891  1.00  4.80           C  
ATOM     60  N   ARG A  10       7.647   4.909  10.005  1.00  3.73           N  
ATOM     61  CA  ARG A  10       8.496   4.609   8.837  1.00  3.38           C  
ATOM     62  C   ARG A  10       7.798   3.609   7.876  1.00  3.47           C  
ATOM     63  O   ARG A  10       7.878   3.778   6.651  1.00  4.67           O  
ATOM     64  CB  ARG A  10       9.847   4.020   9.305  1.00  3.95           C  
ATOM     65  CG  ARG A  10      10.752   3.607   8.149  1.00  4.55           C  
ATOM     66  CD  ARG A  10      11.226   4.699   7.244  1.00  5.89           C  
ATOM     67  NE  ARG A  10      12.143   5.571   8.035  1.00  6.20           N  
ATOM     68  CZ  ARG A  10      12.758   6.609   7.443  1.00  7.52           C  
ATOM     69  NH1 ARG A  10      12.539   6.932   6.158  1.00 10.68           N  
ATOM     70  NH2 ARG A  10      13.601   7.322   8.202  1.00  9.48           N  
ATOM     71  N   SER A  11       7.186   2.582   8.445  1.00  5.19           N  
ATOM     72  CA  SER A  11       6.500   1.584   7.565  1.00  4.60           C  
ATOM     73  C   SER A  11       5.382   2.313   6.773  1.00  4.84           C  
ATOM     74  O   SER A  11       5.213   2.016   5.557  1.00  5.84           O  
ATOM     75  CB  SER A  11       5.908   0.462   8.400  1.00  5.91           C  
ATOM     76  OG  SER A  11       6.990  -0.272   9.012  1.00  8.38           O  
ATOM     77  N   ASN A  12       4.648   3.182   7.446  1.00  3.54           N  
ATOM     78  CA  ASN A  12       3.545   3.935   6.751  1.00  4.57           C  
ATOM     79  C   ASN A  12       4.107   4.851   5.691  1.00  4.14           C  
ATOM     80  O   ASN A  12       3.536   5.001   4.617  1.00  5.52           O  
ATOM     81  CB  ASN A  12       2.663   4.677   7.748  1.00  6.42           C  
ATOM     82  CG  ASN A  12       1.802   3.735   8.610  1.00  8.25           C  
ATOM     83  OD1 ASN A  12       1.567   2.613   8.165  1.00 12.72           O  
ATOM     84  ND2 ASN A  12       1.394   4.252   9.767  1.00  9.92           N  
ATOM     85  N   PHE A  13       5.259   5.498   6.005  1.00  3.43           N  
ATOM     86  CA  PHE A  13       5.929   6.358   5.055  1.00  3.49           C  
ATOM     87  C   PHE A  13       6.304   5.578   3.799  1.00  3.40           C  
ATOM     88  O   PHE A  13       6.136   6.072   2.653  1.00  4.07           O  
ATOM     89  CB  PHE A  13       7.183   6.994   5.754  1.00  5.48           C  
ATOM     90  CG  PHE A  13       7.884   8.006   4.883  1.00  5.57           C  
ATOM     91  CD1 PHE A  13       8.906   7.586   4.027  1.00  6.99           C  
ATOM     92  CD2 PHE A  13       7.532   9.373   4.983  1.00  6.52           C  
ATOM     93  CE1 PHE A  13       9.560   8.539   3.194  1.00  8.20           C  
ATOM     94  CE2 PHE A  13       8.176  10.281   4.145  1.00  6.34           C  
ATOM     95  CZ  PHE A  13       9.141   9.845   3.292  1.00  6.84           C  
ATOM     96  N   ASN A  14       6.900   4.390   3.989  1.00  3.64           N  
ATOM     97  CA  ASN A  14       7.331   3.607   2.791  1.00  4.31           C  
ATOM     98  C   ASN A  14       6.116   3.210   1.915  1.00  3.98           C  
ATOM     99  O   ASN A  14       6.240   3.144   0.684  1.00  6.22           O  
ATOM    100  CB  ASN A  14       8.145   2.404   3.240  1.00  5.81           C  
ATOM    101  CG  ASN A  14       9.555   2.856   3.730  1.00  6.82           C  
ATOM    102  OD1 ASN A  14      10.013   3.895   3.323  1.00  9.43           O  
ATOM    103  ND2 ASN A  14      10.120   1.956   4.539  1.00  8.21           N  
ATOM    104  N   VAL A  15       4.993   2.927   2.571  1.00  3.76           N  
ATOM    105  CA  VAL A  15       3.782   2.599   1.742  1.00  3.98           C  
ATOM    106  C   VAL A  15       3.296   3.871   1.004  1.00  3.80           C  
ATOM    107  O   VAL A  15       2.947   3.817  -0.189  1.00  4.85           O  
ATOM    108  CB  VAL A  15       2.698   1.953   2.608  1.00  4.71           C  
ATOM    109  CG1 VAL A  15       1.384   1.826   1.806  1.00  6.67           C  
ATOM    110  CG2 VAL A  15       3.174   0.533   3.005  1.00  6.26           C  
ATOM    111  N   CYS A  16       3.321   4.987   1.720  1.00  3.79           N  
ATOM    112  CA  CYS A  16       2.890   6.285   1.126  1.00  3.54           C  
ATOM    113  C   CYS A  16       3.687   6.597  -0.111  1.00  3.48           C  
ATOM    114  O   CYS A  16       3.200   7.147  -1.103  1.00  4.63           O  
ATOM    115  CB  CYS A  16       3.039   7.369   2.240  1.00  4.58           C  
ATOM    116  SG  CYS A  16       2.559   9.014   1.649  1.00  5.66           S  
ATOM    117  N   ARG A  17       4.997   6.227  -0.100  1.00  3.99           N  
ATOM    118  CA  ARG A  17       5.895   6.489  -1.213  1.00  3.83           C  
ATOM    119  C   ARG A  17       5.738   5.560  -2.409  1.00  3.79           C  
ATOM    120  O   ARG A  17       6.228   5.901  -3.507  1.00  5.39           O  
ATOM    121  CB  ARG A  17       7.370   6.507  -0.731  1.00  4.11           C  
ATOM    122  CG  ARG A  17       7.717   7.687   0.206  1.00  4.69           C  
ATOM    123  CD  ARG A  17       7.949   8.947  -0.615  1.00  5.10           C  
ATOM    124  NE  ARG A  17       9.212   8.856  -1.337  1.00  4.71           N  
ATOM    125  CZ  ARG A  17       9.537   9.533  -2.431  1.00  5.28           C  
ATOM    126  NH1 ARG A  17       8.659  10.350  -3.032  1.00  6.67           N  
ATOM    127  NH2 ARG A  17      10.793   9.491  -2.899  1.00  6.41           N  
ATOM    128  N   LEU A  18       5.051   4.411  -2.204  1.00  4.70           N  
ATOM    129  CA  LEU A  18       4.933   3.431  -3.326  1.00  5.46           C  
ATOM    130  C   LEU A  18       4.397   4.014  -4.620  1.00  5.13           C  
ATOM    131  O   LEU A  18       4.988   3.755  -5.687  1.00  5.55           O  
ATOM    132  CB  LEU A  18       4.196   2.184  -2.863  1.00  6.47           C  
ATOM    133  CG  LEU A  18       4.960   1.178  -1.991  1.00  7.43           C  
ATOM    134  CD1 LEU A  18       3.907   0.097  -1.634  1.00  8.70           C  
ATOM    135  CD2 LEU A  18       6.129   0.606  -2.768  1.00  9.39           C  
ATOM    136  N   PRO A  19       3.329   4.795  -4.543  1.00  4.28           N  
ATOM    137  CA  PRO A  19       2.792   5.376  -5.797  1.00  5.38           C  
ATOM    138  C   PRO A  19       3.573   6.540  -6.322  1.00  6.30           C  
ATOM    139  O   PRO A  19       3.260   7.045  -7.422  1.00  9.62           O  
ATOM    140  CB  PRO A  19       1.358   5.766  -5.472  1.00  5.87           C  
ATOM    141  CG  PRO A  19       1.223   5.694  -3.993  1.00  6.47           C  
ATOM    142  CD  PRO A  19       2.421   4.941  -3.408  1.00  6.45           C  
ATOM    143  N   GLY A  20       4.565   7.047  -5.559  1.00  4.94           N  
ATOM    144  CA  GLY A  20       5.366   8.191  -6.018  1.00  5.39           C  
ATOM    145  C   GLY A  20       5.007   9.481  -5.280  1.00  5.03           C  
ATOM    146  O   GLY A  20       5.535  10.510  -5.730  1.00  7.34           O  
ATOM    147  N   THR A  21       4.181   9.438  -4.262  1.00  4.10           N  
ATOM    148  CA  THR A  21       3.767  10.609  -3.513  1.00  3.94           C  
ATOM    149  C   THR A  21       5.017  11.397  -3.042  1.00  3.96           C  
ATOM    150  O   THR A  21       5.947  10.757  -2.523  1.00  5.82           O  
ATOM    151  CB  THR A  21       2.992  10.188  -2.225  1.00  4.13           C  
ATOM    152  OG1 THR A  21       2.051   9.144  -2.623  1.00  5.45           O  
ATOM    153  CG2 THR A  21       2.260  11.349  -1.551  1.00  5.41           C  
ATOM    154  N   PRO A  22       4.971  12.703  -3.176  1.00  5.04           N  
ATOM    155  CA  PRO A  22       6.143  13.513  -2.696  1.00  4.69           C  
ATOM    156  C   PRO A  22       6.400  13.233  -1.225  1.00  4.19           C  
ATOM    157  O   PRO A  22       5.485  13.061  -0.382  1.00  4.47           O  
ATOM    158  CB  PRO A  22       5.703  14.969  -2.920  1.00  7.12           C  
ATOM    159  CG  PRO A  22       4.676  14.893  -3.996  1.00  7.03           C  
ATOM    160  CD  PRO A  22       3.964  13.567  -3.811  1.00  4.90           C  
ATOM    161  N   GLU A  23       7.728  13.297  -0.921  1.00  5.16           N  
ATOM    162  CA  GLU A  23       8.114  13.103   0.500  1.00  5.31           C  
ATOM    163  C   GLU A  23       7.427  14.073   1.410  1.00  4.11           C  
ATOM    164  O   GLU A  23       7.036  13.682   2.540  1.00  5.11           O  
ATOM    165  CB  GLU A  23       9.648  13.285   0.660  1.00  6.16           C  
ATOM    166  CG  GLU A  23      10.440  12.093   0.063  1.00  7.48           C  
ATOM    167  CD  GLU A  23      11.941  12.170   0.391  1.00  9.40           C  
ATOM    168  OE1 GLU A  23      12.416  13.225   0.681  1.00 10.40           O  
ATOM    169  OE2 GLU A  23      12.539  11.070   0.292  1.00 13.32           O  
ATOM    170  N   ALA A  24       7.212  15.334   0.966  1.00  4.56           N  
ATOM    171  CA  ALA A  24       6.614  16.317   1.913  1.00  4.49           C  
ATOM    172  C   ALA A  24       5.212  15.936   2.350  1.00  4.10           C  
ATOM    173  O   ALA A  24       4.782  16.166   3.495  1.00  5.64           O  
ATOM    174  CB  ALA A  24       6.605  17.695   1.246  1.00  5.80           C  
ATOM    175  N   ILE A  25       4.445  15.318   1.405  1.00  4.37           N  
ATOM    176  CA  ILE A  25       3.074  14.894   1.756  1.00  5.44           C  
ATOM    177  C   ILE A  25       3.085  13.643   2.645  1.00  4.32           C  
ATOM    178  O   ILE A  25       2.315  13.523   3.578  1.00  4.72           O  
ATOM    179  CB  ILE A  25       2.204  14.637   0.462  1.00  6.42           C  
ATOM    180  CG1 ILE A  25       1.815  16.048  -0.129  1.00  7.50           C  
ATOM    181  CG2 ILE A  25       0.903  13.864   0.811  1.00  7.65           C  
ATOM    182  CD1 ILE A  25       0.756  16.761   0.757  1.00  7.80           C  
ATOM    183  N   CYS A  26       4.032  12.764   2.313  1.00  3.92           N  
ATOM    184  CA  CYS A  26       4.180  11.549   3.187  1.00  4.37           C  
ATOM    185  C   CYS A  26       4.632  11.944   4.596  1.00  3.95           C  
ATOM    186  O   CYS A  26       4.227  11.252   5.547  1.00  4.74           O  
ATOM    187  CB  CYS A  26       5.038  10.518   2.539  1.00  4.63           C  
ATOM    188  SG  CYS A  26       4.349   9.794   1.022  1.00  5.61           S  
ATOM    189  N   ALA A  27       5.408  13.012   4.694  1.00  3.89           N  
ATOM    190  CA  ALA A  27       5.879  13.502   6.026  1.00  4.43           C  
ATOM    191  C   ALA A  27       4.696  13.908   6.882  1.00  4.26           C  
ATOM    192  O   ALA A  27       4.528  13.422   8.025  1.00  5.44           O  
ATOM    193  CB  ALA A  27       6.880  14.615   5.830  1.00  5.36           C  
ATOM    194  N   THR A  28       3.827  14.802   6.358  1.00  4.53           N  
ATOM    195  CA  THR A  28       2.691  15.221   7.194  1.00  5.08           C  
ATOM    196  C   THR A  28       1.672  14.132   7.434  1.00  4.62           C  
ATOM    197  O   THR A  28       0.947  14.112   8.468  1.00  7.80           O  
ATOM    198  CB  THR A  28       1.986  16.520   6.614  1.00  6.03           C  
ATOM    199  OG1 THR A  28       1.664  16.221   5.230  1.00  7.19           O  
ATOM    200  CG2 THR A  28       2.914  17.739   6.700  1.00  7.34           C  
ATOM    201  N   TYR A  29       1.621  13.190   6.511  1.00  5.01           N  
ATOM    202  CA  TYR A  29       0.715  12.045   6.657  1.00  6.60           C  
ATOM    203  C   TYR A  29       1.125  11.125   7.815  1.00  4.92           C  
ATOM    204  O   TYR A  29       0.286  10.632   8.545  1.00  7.13           O  
ATOM    205  CB  TYR A  29       0.755  11.229   5.322  1.00  9.66           C  
ATOM    206  CG  TYR A  29      -0.203  10.044   5.354  1.00 11.56           C  
ATOM    207  CD1 TYR A  29      -1.547  10.337   5.645  1.00 12.85           C  
ATOM    208  CD2 TYR A  29       0.193   8.750   5.100  1.00 14.44           C  
ATOM    209  CE1 TYR A  29      -2.496   9.329   5.673  1.00 16.61           C  
ATOM    210  CE2 TYR A  29      -0.801   7.705   5.156  1.00 17.11           C  
ATOM    211  CZ  TYR A  29      -2.079   8.031   5.430  1.00 19.99           C  
ATOM    212  OH  TYR A  29      -3.097   7.057   5.458  1.00 28.98           O  
ATOM    213  N   THR A  30       2.470  10.984   7.995  1.00  5.31           N  
ATOM    214  CA  THR A  30       2.986   9.994   8.950  1.00  5.70           C  
ATOM    215  C   THR A  30       3.609  10.505  10.230  1.00  6.28           C  
ATOM    216  O   THR A  30       3.766   9.715  11.186  1.00  8.77           O  
ATOM    217  CB  THR A  30       4.076   9.103   8.225  1.00  6.55           C  
ATOM    218  OG1 THR A  30       5.125  10.027   7.824  1.00  6.57           O  
ATOM    219  CG2 THR A  30       3.493   8.324   7.035  1.00  7.29           C  
ATOM    220  N   GLY A  31       3.984  11.764  10.241  1.00  4.99           N  
ATOM    221  CA  GLY A  31       4.769  12.336  11.360  1.00  5.50           C  
ATOM    222  C   GLY A  31       6.255  12.243  11.106  1.00  4.19           C  
ATOM    223  O   GLY A  31       7.037  12.750  11.954  1.00  6.12           O  
ATOM    224  N   CYS A  32       6.710  11.631   9.992  1.00  4.30           N  
ATOM    225  CA  CYS A  32       8.140  11.694   9.635  1.00  4.89           C  
ATOM    226  C   CYS A  32       8.500  13.141   9.206  1.00  5.50           C  
ATOM    227  O   CYS A  32       7.581  13.949   8.944  1.00  5.82           O  
ATOM    228  CB  CYS A  32       8.504  10.686   8.530  1.00  4.66           C  
ATOM    229  SG  CYS A  32       8.048   8.987   8.881  1.00  5.33           S  
ATOM    230  N   ILE A  33       9.793  13.410   9.173  1.00  6.02           N  
ATOM    231  CA  ILE A  33      10.280  14.760   8.823  1.00  5.24           C  
ATOM    232  C   ILE A  33      11.346  14.658   7.743  1.00  5.16           C  
ATOM    233  O   ILE A  33      11.971  13.583   7.552  1.00  7.19           O  
ATOM    234  CB  ILE A  33      10.790  15.535  10.085  1.00  5.49           C  
ATOM    235  CG1 ILE A  33      12.059  14.803  10.671  1.00  6.85           C  
ATOM    236  CG2 ILE A  33       9.684  15.686  11.138  1.00  6.45           C  
ATOM    237  CD1 ILE A  33      12.733  15.676  11.781  1.00  8.94           C  
ATOM    238  N   ILE A  34      11.490  15.773   7.038  1.00  5.52           N  
ATOM    239  CA  ILE A  34      12.552  15.877   6.036  1.00  6.82           C  
ATOM    240  C   ILE A  34      13.590  16.917   6.560  1.00  6.92           C  
ATOM    241  O   ILE A  34      13.168  18.006   6.945  1.00  9.22           O  
ATOM    242  CB  ILE A  34      11.987  16.360   4.681  1.00  8.11           C  
ATOM    243  CG1 ILE A  34      10.914  15.338   4.163  1.00  9.59           C  
ATOM    244  CG2 ILE A  34      13.131  16.517   3.629  1.00  9.73           C  
ATOM    245  CD1 ILE A  34      10.151  16.024   2.938  1.00 13.41           C  
ATOM    246  N   ILE A  35      14.856  16.493   6.536  1.00  7.06           N  
ATOM    247  CA  ILE A  35      15.930  17.454   6.941  1.00  7.52           C  
ATOM    248  C   ILE A  35      16.913  17.550   5.819  1.00  6.63           C  
ATOM    249  O   ILE A  35      17.097  16.660   4.970  1.00  7.90           O  
ATOM    250  CB  ILE A  35      16.622  16.995   8.285  1.00  8.07           C  
ATOM    251  CG1 ILE A  35      17.360  15.651   8.067  1.00  9.41           C  
ATOM    252  CG2 ILE A  35      15.592  16.974   9.434  1.00  9.46           C  
ATOM    253  CD1 ILE A  35      18.298  15.206   9.219  1.00  9.85           C  
ATOM    254  N   PRO A  36      17.664  18.669   5.806  1.00  8.07           N  
ATOM    255  CA  PRO A  36      18.635  18.861   4.738  1.00  8.78           C  
ATOM    256  C   PRO A  36      19.925  18.042   4.949  1.00  8.31           C  
ATOM    257  O   PRO A  36      20.593  17.742   3.945  1.00  9.09           O  
ATOM    258  CB  PRO A  36      18.945  20.364   4.783  1.00  9.67           C  
ATOM    259  CG  PRO A  36      18.238  20.937   5.908  1.00 10.15           C  
ATOM    260  CD  PRO A  36      17.371  19.900   6.596  1.00  9.53           C  
ATOM    261  N   GLY A  37      20.172  17.730   6.217  1.00  8.48           N  
ATOM    262  CA  GLY A  37      21.452  16.969   6.513  1.00  9.20           C  
ATOM    263  C   GLY A  37      21.143  15.478   6.427  1.00 10.41           C  
ATOM    264  O   GLY A  37      20.138  15.023   5.878  1.00 12.06           O  
ATOM    265  N   ALA A  38      22.055  14.701   7.032  1.00  9.24           N  
ATOM    266  CA  ALA A  38      22.019  13.242   7.020  1.00  9.24           C  
ATOM    267  C   ALA A  38      21.944  12.628   8.396  1.00  9.60           C  
ATOM    268  O   ALA A  38      21.869  11.387   8.435  1.00 13.65           O  
ATOM    269  CB  ALA A  38      23.246  12.697   6.275  1.00 10.43           C  
ATOM    270  N   THR A  39      21.894  13.435   9.436  1.00  8.70           N  
ATOM    271  CA  THR A  39      21.936  12.911  10.809  1.00  9.46           C  
ATOM    272  C   THR A  39      20.615  13.191  11.521  1.00  8.32           C  
ATOM    273  O   THR A  39      20.357  14.317  11.948  1.00  9.89           O  
ATOM    274  CB  THR A  39      23.131  13.601  11.593  1.00 10.72           C  
ATOM    275  OG1 THR A  39      24.284  13.401  10.709  1.00 11.66           O  
ATOM    276  CG2 THR A  39      23.340  12.935  12.962  1.00 11.81           C  
ATOM    277  N   CYS A  40      19.827  12.110  11.642  1.00  7.64           N  
ATOM    278  CA  CYS A  40      18.504  12.312  12.298  1.00  8.05           C  
ATOM    279  C   CYS A  40      18.684  12.451  13.784  1.00  7.63           C  
ATOM    280  O   CYS A  40      19.533  11.718  14.362  1.00  9.64           O  
ATOM    281  CB  CYS A  40      17.582  11.117  11.996  1.00  7.80           C  
ATOM    282  SG  CYS A  40      17.199  10.929  10.237  1.00  7.30           S  
ATOM    283  N   PRO A  41      17.880  13.266  14.426  1.00  8.00           N  
ATOM    284  CA  PRO A  41      17.924  13.421  15.877  1.00  8.96           C  
ATOM    285  C   PRO A  41      17.392  12.206  16.594  1.00  9.06           C  
ATOM    286  O   PRO A  41      16.652  11.368  16.033  1.00  8.82           O  
ATOM    287  CB  PRO A  41      17.076  14.658  16.145  1.00 10.39           C  
ATOM    288  CG  PRO A  41      16.098  14.689  14.997  1.00 10.99           C  
ATOM    289  CD  PRO A  41      16.859  14.150  13.779  1.00 10.49           C  
ATOM    290  N   GLY A  42      17.728  12.124  17.884  1.00  7.55           N  
ATOM    291  CA  GLY A  42      17.334  10.956  18.691  1.00  8.00           C  
ATOM    292  C   GLY A  42      15.875  10.688  18.871  1.00  7.22           C  
ATOM    293  O   GLY A  42      15.434   9.550  19.166  1.00  8.41           O  
ATOM    294  N   ASP A  43      15.036  11.747  18.715  1.00  5.54           N  
ATOM    295  CA  ASP A  43      13.564  11.573  18.836  1.00  5.85           C  
ATOM    296  C   ASP A  43      12.936  11.227  17.470  1.00  5.87           C  
ATOM    297  O   ASP A  43      11.720  11.040  17.428  1.00  7.29           O  
ATOM    298  CB  ASP A  43      12.933  12.737  19.580  1.00  6.72           C  
ATOM    299  CG  ASP A  43      13.140  14.094  18.958  1.00  8.59           C  
ATOM    300  OD1 ASP A  43      14.109  14.303  18.212  1.00  9.59           O  
ATOM    301  OD2 ASP A  43      12.267  14.963  19.265  1.00 11.45           O  
ATOM    302  N   TYR A  44      13.725  11.174  16.425  1.00  5.22           N  
ATOM    303  CA  TYR A  44      13.257  10.745  15.081  1.00  5.56           C  
ATOM    304  C   TYR A  44      14.275   9.687  14.612  1.00  4.61           C  
ATOM    305  O   TYR A  44      14.930   9.862  13.568  1.00  6.04           O  
ATOM    306  CB  TYR A  44      13.200  11.914  14.071  1.00  5.41           C  
ATOM    307  CG  TYR A  44      12.000  12.819  14.399  1.00  5.34           C  
ATOM    308  CD1 TYR A  44      12.119  13.853  15.332  1.00  6.59           C  
ATOM    309  CD2 TYR A  44      10.775  12.617  13.762  1.00  5.94           C  
ATOM    310  CE1 TYR A  44      11.045  14.675  15.610  1.00  5.97           C  
ATOM    311  CE2 TYR A  44       9.676  13.433  14.048  1.00  5.17           C  
ATOM    312  CZ  TYR A  44       9.802  14.456  14.996  1.00  5.96           C  
ATOM    313  OH  TYR A  44       8.740  15.265  15.269  1.00  8.60           O  
ATOM    314  N   ALA A  45      14.342   8.640  15.422  1.00  4.76           N  
ATOM    315  CA  ALA A  45      15.445   7.667  15.246  1.00  5.89           C  
ATOM    316  C   ALA A  45      15.171   6.533  14.280  1.00  6.67           C  
ATOM    317  O   ALA A  45      16.093   5.705  14.039  1.00  7.56           O  
ATOM    318  CB  ALA A  45      15.680   7.099  16.682  1.00  6.82           C  
ATOM    319  N   ASN A  46      13.966   6.502  13.739  1.00  5.80           N  
ATOM    320  CA  ASN A  46      13.512   5.395  12.878  1.00  6.15           C  
ATOM    321  C   ASN A  46      13.311   5.853  11.455  1.00  6.61           C  
ATOM    322  O   ASN A  46      13.733   6.929  11.026  1.00  7.18           O  
ATOM    323  CB  ASN A  46      12.266   4.769  13.501  1.00  7.27           C  
ATOM    324  CG  ASN A  46      12.538   4.304  14.922  1.00  7.98           C  
ATOM    325  OD1 ASN A  46      11.982   4.849  15.886  1.00 11.00           O  
ATOM    326  ND2 ASN A  46      13.407   3.298  15.015  1.00 10.32           N  
ATOM    327  OXT ASN A  46      12.703   4.973  10.746  1.00  7.86           O  
CONECT    1    2
CONECT    2    3    5
CONECT    3    4    4    8
CONECT    5    6    7
CONECT    8    9
CONECT    9   10   12
CONECT   10   11   11   15
CONECT   12   13   14
CONECT   15   16
CONECT   16   17   19
CONECT   17   18   18   21
CONECT   19   20
CONECT   20  282
CONECT   21   22
CONECT   22   23   25
CONECT   23   24   24   27
CONECT   25   26
CONECT   26  229
CONECT   27   28   33
CONECT   28   29   31
CONECT   29   30   30   34
CONECT   31   32
CONECT   32   33
CONECT   34   35
CONECT   35   36   38
CONECT   36   37   37   40
CONECT   38   39
CONECT   40   41
CONECT   41   42   44
CONECT   42   43   43   48
CONECT   44   45   46
CONECT   45   47
CONECT   48   49
CONECT   49   50   52
CONECT   50   51   51   55
CONECT   52   53   54
CONECT   55   56
CONECT   56   57   59
CONECT   57   58   58   60
CONECT   60   61
CONECT   61   62   64
CONECT   62   63   63   71
CONECT   64   65
CONECT   65   66
CONECT   66   67
CONECT   67   68
CONECT   68   69   70   70
CONECT   71   72
CONECT   72   73   75
CONECT   73   74   74   77
CONECT   75   76
CONECT   77   78
CONECT   78   79   81
CONECT   79   80   80   85
CONECT   81   82
CONECT   82   83   83   84
CONECT   85   86
CONECT   86   87   89
CONECT   87   88   88   96
CONECT   89   90
CONECT   90   91   91   92
CONECT   91   93
CONECT   92   94   94
CONECT   93   95   95
CONECT   94   95
CONECT   96   97
CONECT   97   98  100
CONECT   98   99   99  104
CONECT  100  101
CONECT  101  102  102  103
CONECT  104  105
CONECT  105  106  108
CONECT  106  107  107  111
CONECT  108  109  110
CONECT  111  112
CONECT  112  113  115
CONECT  113  114  114  117
CONECT  115  116
CONECT  116  188
CONECT  117  118
CONECT  118  119  121
CONECT  119  120  120  128
CONECT  121  122
CONECT  122  123
CONECT  123  124
CONECT  124  125
CONECT  125  126  127  127
CONECT  128  129
CONECT  129  130  132
CONECT  130  131  131  136
CONECT  132  133
CONECT  133  134  135
CONECT  136  137  142
CONECT  137  138  140
CONECT  138  139  139  143
CONECT  140  141
CONECT  141  142
CONECT  143  144
CONECT  144  145
CONECT  145  146  146  147
CONECT  147  148
CONECT  148  149  151
CONECT  149  150  150  154
CONECT  151  152  153
CONECT  154  155  160
CONECT  155  156  158
CONECT  156  157  157  161
CONECT  158  159
CONECT  159  160
CONECT  161  162
CONECT  162  163  165
CONECT  163  164  164  170
CONECT  165  166
CONECT  166  167
CONECT  167  168  168  169
CONECT  170  171
CONECT  171  172  174
CONECT  172  173  173  175
CONECT  175  176
CONECT  176  177  179
CONECT  177  178  178  183
CONECT  179  180  181
CONECT  180  182
CONECT  183  184
CONECT  184  185  187
CONECT  185  186  186  189
CONECT  187  188
CONECT  189  190
CONECT  190  191  193
CONECT  191  192  192  194
CONECT  194  195
CONECT  195  196  198
CONECT  196  197  197  201
CONECT  198  199  200
CONECT  201  202
CONECT  202  203  205
CONECT  203  204  204  213
CONECT  205  206
CONECT  206  207  207  208
CONECT  207  209
CONECT  208  210  210
CONECT  209  211  211
CONECT  210  211
CONECT  211  212
CONECT  213  214
CONECT  214  215  217
CONECT  215  216  216  220
CONECT  217  218  219
CONECT  220  221
CONECT  221  222
CONECT  222  223  223  224
CONECT  224  225
CONECT  225  226  228
CONECT  226  227  227  230
CONECT  228  229
CONECT  230  231
CONECT  231  232  234
CONECT  232  233  233  238
CONECT  234  235  236
CONECT  235  237
CONECT  238  239
CONECT  239  240  242
CONECT  240  241  241  246
CONECT  242  243  244
CONECT  243  245
CONECT  246  247
CONECT  247  248  250
CONECT  248  249  249  254
CONECT  250  251  252
CONECT  251  253
CONECT  254  255  260
CONECT  255  256  258
CONECT  256  257  257  261
CONECT  258  259
CONECT  259  260
CONECT  261  262
CONECT  262  263
CONECT  263  264  264  265
CONECT  265  266
CONECT  266  267  269
CONECT  267  268  268  270
CONECT  270  271
CONECT  271  272  274
CONECT  272  273  273  277
CONECT  274  275  276
CONECT  277  278
CONECT  278  279  281
CONECT  279  280  280  283
CONECT  281  282
CONECT  283  284  289
CONECT  284  285  287
CONECT  285  286  286  290
CONECT  287  288
CONECT  288  289
CONECT  290  291
CONECT  291  292
CONECT  292  293  293  294
CONECT  294  295
CONECT  295  296  298
CONECT  296  297  297  302
CONECT  298  299
CONECT  299  300  300  301
CONECT  302  303
CONECT  303  304  306
CONECT  304  305  305  314
CONECT  306  307
CONECT  307  308  308  309
CONECT  308  310
CONECT  309  311  311
CONECT  310  312  312
CONECT  311  312
CONECT  312  313
CONECT  314  315
CONECT  315  316  318
CONECT  316  317  317  319
CONECT  319  320
CONECT  320  321  323
CONECT  321  322  322  327
CONECT  323  324
CONECT  324  325  325  326
END
</textarea>*/
