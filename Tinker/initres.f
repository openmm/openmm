c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine initres  --  setup biopolymer residue names  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "initres" sets biopolymer residue names and biotype codes used
c     in PDB file conversion and automated generation of structures
c
c
      subroutine initres
      use resdue
      implicit none
      integer i
      integer nt(maxamino),cat(maxamino)
      integer ct(maxamino),hnt(maxamino)
      integer ot(maxamino),hat(maxamino)
      integer cbt(maxamino)
      integer nn(maxamino),can(maxamino)
      integer cn(maxamino),hnn(maxamino)
      integer on(maxamino),han(maxamino)
      integer nc(maxamino),cac(maxamino)
      integer cc(maxamino),hnc(maxamino)
      integer oc(maxamino),hac(maxamino)
      integer o5t(maxnuc),c5t(maxnuc)
      integer h51t(maxnuc),h52t(maxnuc)
      integer c4t(maxnuc),h4t(maxnuc)
      integer o4t(maxnuc),c1t(maxnuc)
      integer h1t(maxnuc),c3t(maxnuc)
      integer h3t(maxnuc),c2t(maxnuc)
      integer h21t(maxnuc),o2t(maxnuc)
      integer h22t(maxnuc),o3t(maxnuc)
      integer pt(maxnuc),opt(maxnuc)
      integer h5tt(maxnuc),h3tt(maxnuc)
      character*1 acid1(maxamino)
      character*1 base1(maxnuc)
      character*3 acid3(maxamino)
      character*3 base3(maxnuc)
c
c     supported amino acid 1-letter and 3-letter codes
c
      data acid1  / 'G', 'A', 'V', 'L', 'I', 'S', 'T', 'C', 'C', 'c',
     &              'P', 'F', 'Y', 'y', 'W', 'H', 'U', 'Z', 'D', 'd',
     &              'N', 'E', 'e', 'Q', 'M', 'K', 'k', 'R', 'O', 'B',
     &              'J', 't', 'f', 'a', 'x', 'n', 'm', 'X' /
      data acid3  / 'GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'SER', 'THR',
     &              'CYS', 'CYX', 'CYD', 'PRO', 'PHE', 'TYR', 'TYD',
     &              'TRP', 'HIS', 'HID', 'HIE', 'ASP', 'ASH', 'ASN',
     &              'GLU', 'GLH', 'GLN', 'MET', 'LYS', 'LYD', 'ARG',
     &              'ORN', 'AIB', 'PCA', 'H2N', 'FOR', 'ACE', 'COH',
     &              'NH2', 'NME', 'UNK' /
c
c     supported nucleotide 1-letter and 3-letter codes
c
      data base1  / 'A', 'G', 'C', 'U', 'D', 'B', 'I', 'T', '1', '2',
     &              '3', 'X' /
      data base3  / '  A', '  G', '  C', '  U', ' DA', ' DG', ' DC',
     &              ' DT', ' MP', ' DP', ' TP', 'UNK' /
c
c     biopolymer types for mid-chain peptide backbone atoms
c
      data nt   /   1,   7,  15,  27,  41,  55,  65,  77,  87,  96,
     &            105, 116, 131, 147, 162, 185, 202, 218, 234, 244,
     &            256, 268, 280, 294, 308, 321, 337, 353, 370, 384,
     &            391,   0,   0,   0,   0,   0,   0,   1 /
      data cat  /   2,   8,  16,  28,  42,  56,  66,  78,  88,  97,
     &            106, 117, 132, 148, 163, 186, 203, 219, 235, 245,
     &            257, 269, 281, 295, 309, 322, 338, 354, 371, 385,
     &            392,   0,   0,   0,   0,   0,   0,   2 /
      data ct   /   3,   9,  17,  29,  43,  57,  67,  79,  89,  98,
     &            107, 118, 133, 149, 164, 187, 204, 220, 236, 246,
     &            258, 270, 282, 296, 310, 323, 339, 355, 372, 386,
     &            393,   0,   0,   0,   0,   0,   0,   3 /
      data hnt  /   4,  10,  18,  30,  44,  58,  68,  80,  90,  99,
     &              0, 119, 134, 150, 165, 188, 205, 221, 237, 247,
     &            259, 271, 283, 297, 311, 324, 340, 356, 373, 387,
     &            394,   0,   0,   0,   0,   0,   0,   4 /
      data ot   /   5,  11,  19,  31,  45,  59,  69,  81,  91, 100,
     &            108, 120, 135, 151, 166, 189, 206, 222, 238, 248,
     &            260, 272, 284, 298, 312, 325, 341, 357, 374, 388,
     &            395,   0,   0,   0,   0,   0,   0,   5 /
      data hat  /   6,  12,  20,  32,  46,  60,  70,  82,  92, 101,
     &            109, 121, 136, 152, 167, 190, 207, 223, 239, 249,
     &            261, 273, 285, 299, 313, 326, 342, 358, 375,   0,
     &            396,   0,   0,   0,   0,   0,   0,   6 /
      data cbt  /   0,  13,  21,  33,  47,  61,  71,  83,  93, 102,
     &            110, 122, 137, 153, 168, 191, 208, 224, 240, 250,
     &            262, 274, 286, 300, 314, 327, 343, 359, 376, 389,
     &            397,   0,   0,   0,   0,   0,   0,   0 /
c
c     biopolymer types for N-terminal peptide backbone atoms
c
      data nn   / 403, 409, 415, 421, 427, 433, 439, 445, 451, 457,
     &            463, 471, 477, 483, 489, 495, 501, 507, 513, 519,
     &            525, 531, 537, 543, 549, 555, 561, 567, 573, 579,
     &            391, 762,   0,   0,   0,   0,   0, 403 /
      data can  / 404, 410, 416, 422, 428, 434, 440, 446, 452, 458,
     &            464, 472, 478, 484, 490, 496, 502, 508, 514, 520,
     &            526, 532, 538, 544, 550, 556, 562, 568, 574, 580,
     &            392,   0,   0, 767,   0,   0,   0, 404 /
      data cn   / 405, 411, 417, 423, 429, 435, 441, 447, 453, 459,
     &            465, 473, 479, 485, 491, 497, 503, 509, 515, 521,
     &            527, 533, 539, 545, 551, 557, 563, 569, 575, 581,
     &            393,   0, 764, 769,   0,   0,   0, 405 /
      data hnn  / 406, 412, 418, 424, 430, 436, 442, 448, 454, 460,
     &            466, 474, 480, 486, 492, 498, 504, 510, 516, 522,
     &            528, 534, 540, 546, 552, 558, 564, 570, 576, 582,
     &            394, 763,   0,   0,   0,   0,   0, 406 /
      data on   / 407, 413, 419, 425, 431, 437, 443, 449, 455, 461,
     &            467, 475, 481, 487, 493, 499, 505, 511, 517, 523,
     &            529, 535, 541, 547, 553, 559, 565, 571, 577, 583,
     &            395,   0, 766, 770,   0,   0,   0, 407 /
      data han  / 408, 414, 420, 426, 432, 438, 444, 450, 456, 462,
     &            468, 476, 482, 488, 494, 500, 506, 512, 518, 524,
     &            530, 536, 542, 548, 554, 560, 566, 572, 578,   0,
     &            396,   0, 765, 768,   0,   0,   0, 408 /
c
c     biopolymer types for C-terminal peptide backbone atoms
c
      data nc   / 584, 590, 596, 602, 608, 614, 620, 626, 632, 638,
     &            644, 649, 655, 661, 667, 673, 679, 685, 691, 697,
     &            703, 709, 715, 721, 727, 733, 739, 745, 751, 757,
     &              0,   0,   0,   0, 773, 775, 777, 584 /
      data cac  / 585, 591, 597, 603, 609, 615, 621, 627, 633, 639,
     &            645, 650, 656, 662, 668, 674, 680, 686, 692, 698,
     &            704, 710, 716, 722, 728, 734, 740, 746, 752, 758,
     &              0,   0,   0,   0,   0,   0, 779, 585 /
      data cc   / 586, 592, 598, 604, 610, 616, 622, 628, 634, 640,
     &            646, 651, 657, 663, 669, 675, 681, 687, 693, 699,
     &            705, 711, 717, 723, 729, 735, 741, 747, 753, 759,
     &              0,   0,   0,   0, 771,   0,   0, 586 /
      data hnc  / 587, 593, 599, 605, 611, 617, 623, 629, 635, 641,
     &              0, 652, 658, 664, 670, 676, 682, 688, 694, 700,
     &            706, 712, 718, 724, 730, 736, 742, 748, 754, 760,
     &              0,   0,   0,   0, 774, 776, 778, 587 /
      data oc   / 588, 594, 600, 606, 612, 618, 624, 630, 636, 642,
     &            647, 653, 659, 665, 671, 677, 683, 689, 695, 701,
     &            707, 713, 719, 725, 731, 737, 743, 749, 755, 761,
     &              0,   0,   0,   0, 772,   0,   0, 588 /
      data hac  / 589, 595, 601, 607, 613, 619, 625, 631, 637, 643,
     &            648, 654, 660, 666, 672, 678, 684, 690, 696, 702,
     &            708, 714, 720, 726, 732, 738, 744, 750, 756,   0,
     &              0,   0,   0,   0,   0,   0, 780, 589 /
c
c     biopolymer types for nucleotide phosphate and sugar atoms
c
      data o5t   / 1001, 1031, 1062, 1090, 1117, 1146, 1176, 1203,
     &                0,    0,    0,    0 /
      data c5t   / 1002, 1032, 1063, 1091, 1118, 1147, 1177, 1204,
     &                0,    0,    0,    0 /
      data h51t  / 1003, 1033, 1064, 1092, 1119, 1148, 1178, 1205,
     &                0,    0,    0,    0 /
      data h52t  / 1004, 1034, 1065, 1093, 1120, 1149, 1179, 1206,
     &                0,    0,    0,    0 /
      data c4t   / 1005, 1035, 1066, 1094, 1121, 1150, 1180, 1207,
     &                0,    0,    0,    0 /
      data h4t   / 1006, 1036, 1067, 1095, 1122, 1151, 1181, 1208,
     &                0,    0,    0,    0 /
      data o4t   / 1007, 1037, 1068, 1096, 1123, 1152, 1182, 1209,
     &                0,    0,    0,    0 /
      data c1t   / 1008, 1038, 1069, 1097, 1124, 1153, 1183, 1210,
     &                0,    0,    0,    0 /
      data h1t   / 1009, 1039, 1070, 1098, 1125, 1154, 1184, 1211,
     &                0,    0,    0,    0 /
      data c3t   / 1010, 1040, 1071, 1099, 1126, 1155, 1185, 1212,
     &                0,    0,    0,    0 /
      data h3t   / 1011, 1041, 1072, 1100, 1127, 1156, 1186, 1213,
     &                0,    0,    0,    0 /
      data c2t   / 1012, 1042, 1073, 1101, 1128, 1157, 1187, 1214,
     &                0,    0,    0,    0 /
      data h21t  / 1013, 1043, 1074, 1102, 1129, 1158, 1188, 1215,
     &                0,    0,    0,    0 /
      data o2t   / 1014, 1044, 1075, 1103,    0,    0,    0,    0,
     &                0,    0,    0,    0 /
      data h22t  / 1015, 1045, 1076, 1104, 1130, 1159, 1189, 1216,
     &                0,    0,    0,    0 /
      data o3t   / 1016, 1046, 1077, 1105, 1131, 1160, 1190, 1217,
     &                0,    0,    0,    0 /
      data pt    / 1230, 1230, 1230, 1230, 1242, 1242, 1242, 1242,
     &                0,    0,    0,    0 /
      data opt   / 1231, 1231, 1231, 1231, 1243, 1243, 1243, 1243,
     &                0,    0,    0,    0 /
      data h5tt  / 1233, 1233, 1233, 1233, 1245, 1245, 1245, 1245,
     &                0,    0,    0,    0 /
      data h3tt  / 1238, 1238, 1238, 1238, 1250, 1250, 1250, 1250,
     &                0,    0,    0,    0 /
c
c
c     set amino acid names and peptide backbone biotypes
c
      do i = 1, maxamino
         amino(i) = acid3(i)
         amino1(i) = acid1(i)
         ntyp(i) = nt(i)
         catyp(i) = cat(i)
         ctyp(i) = ct(i)
         hntyp(i) = hnt(i)
         otyp(i) = ot(i)
         hatyp(i) = hat(i)
         cbtyp(i) = cbt(i)
         nntyp(i) = nn(i)
         cantyp(i) = can(i)
         cntyp(i) = cn(i)
         hnntyp(i) = hnn(i)
         ontyp(i) = on(i)
         hantyp(i) = han(i)
         nctyp(i) = nc(i)
         cactyp(i) = cac(i)
         cctyp(i) = cc(i)
         hnctyp(i) = hnc(i)
         octyp(i) = oc(i)
         hactyp(i) = hac(i)
      end do
c
c     set values for the 1- and 3-letter nucleotide names
c
      do i = 1, maxnuc
         nuclz(i) = base3(i)
         nuclz1(i) = base1(i)
         o5typ(i) = o5t(i)
         c5typ(i) = c5t(i)
         h51typ(i) = h51t(i)
         h52typ(i) = h52t(i)
         c4typ(i) = c4t(i)
         h4typ(i) = h4t(i)
         o4typ(i) = o4t(i)
         c1typ(i) = c1t(i)
         h1typ(i) = h1t(i)
         c3typ(i) = c3t(i)
         h3typ(i) = h3t(i)
         c2typ(i) = c2t(i)
         h21typ(i) = h21t(i)
         o2typ(i) = o2t(i)
         h22typ(i) = h22t(i)
         o3typ(i) = o3t(i)
         ptyp(i) = pt(i)
         optyp(i) = opt(i)
         h5ttyp(i) = h5tt(i)
         h3ttyp(i) = h3tt(i)
      end do
      return
      end
