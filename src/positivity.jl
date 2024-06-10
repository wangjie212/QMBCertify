function posepsd8!(model, cons, tsupp, L; lattice="chain")
    ltsupp = length(tsupp)
    Pauli = Matrix{Complex{Int8}}[[1 0; 0 1], [0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]
    blocks = [[1, 16, 24, 28, 30, 31, 40, 44, 46, 47, 52, 54, 55, 58, 59, 61, 72, 76, 78, 79, 84, 86, 87, 90, 91, 93, 100, 102, 103, 106, 107, 109, 114, 115, 117, 121, 136, 140, 142, 143, 148, 150, 151, 154, 155, 157, 164, 166, 167, 170, 171, 173, 178, 179, 181, 185, 196, 198, 199, 202, 203, 205, 210, 211, 213, 217, 226, 227, 229, 233, 241, 256], 
    [2, 3, 5, 9, 17, 32, 33, 48, 56, 60, 62, 63, 65, 80, 88, 92, 94, 95, 104, 108, 110, 111, 116, 118, 119, 122, 123, 125, 129, 144, 152, 156, 158, 159, 168, 172, 174, 175, 180, 182, 183, 186, 187, 189, 200, 204, 206, 207, 212, 214, 215, 218, 219, 221, 228, 230, 231, 234, 235, 237, 242, 243, 245, 249], 
    [4, 6, 7, 10, 11, 13, 18, 19, 21, 25, 34, 35, 37, 41, 49, 64, 66, 67, 69, 73, 81, 96, 97, 112, 120, 124, 126, 127, 130, 131, 133, 137, 145, 160, 161, 176, 184, 188, 190, 191, 193, 208, 216, 220, 222, 223, 232, 236, 238, 239, 244, 246, 247, 250, 251, 253], 
    [8, 12, 14, 15, 20, 22, 23, 26, 27, 29, 36, 38, 39, 42, 43, 45, 50, 51, 53, 57, 68, 70, 71, 74, 75, 77, 82, 83, 85, 89, 98, 99, 101, 105, 113, 128, 132, 134, 135, 138, 139, 141, 146, 147, 149, 153, 162, 163, 165, 169, 177, 192, 194, 195, 197, 201, 209, 224, 225, 240, 248, 252, 254, 255]]
    pos = [@variable(model, [1:length(blocks[i]), 1:length(blocks[i])], PSD) for i =1:4]
    for i = 0:3, j = 0:3, k = 0:3, l = 0:3, s = 0:3, t = 0:3, u = 0:3, v = 0:3
        ind = [i,j,k,l,s,t,u,v]
        if all(x->iseven(sum(ind .== x)), 1:3)
            inx = ind .!= 0
            Locb = bfind(tsupp, ltsupp, reduce4(UInt16.(3*(Vector(1:8)[inx] .- 1) + ind[inx]), L, lattice=lattice))
            tp = real(kron(Pauli[i+1], Pauli[j+1], Pauli[k+1], Pauli[l+1], Pauli[s+1], Pauli[t+1], Pauli[u+1], Pauli[v+1]))
            for p = 1:4
                @inbounds add_to_expression!(cons[Locb], sum(tp[blocks[p], blocks[p]].*pos[p]))
            end
        end
    end
end

function posepsd9!(model, cons, tsupp, L; lattice="chain")
    ltsupp = length(tsupp)
    Pauli = Matrix{Complex{Int8}}[[1 0; 0 1], [0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]
    blocks = [[1, 16, 24, 28, 30, 31, 40, 44, 46, 47, 52, 54, 55, 58, 59, 61, 72, 76, 78, 79, 84, 86, 87, 90, 91, 93, 100, 102, 103, 106, 107, 109, 114, 115, 117, 121, 136, 140, 142, 143, 148, 150, 151, 154, 155, 157, 164, 166, 167, 170, 171, 173, 178, 179, 181, 185, 196, 198, 199, 202, 203, 205, 210, 211, 213, 217, 226, 227, 229, 233, 241, 256, 264, 268, 270, 271, 276, 278, 279, 282, 283, 285, 292, 294, 295, 298, 299, 301, 306, 307, 309, 313, 324, 326, 327, 330, 331, 333, 338, 339, 341, 345, 354, 355, 357, 361, 369, 384, 388, 390, 391, 394, 395, 397, 402, 403, 405, 409, 418, 419, 421, 425, 433, 448, 450, 451, 453, 457, 465, 480, 481, 496, 504, 508, 510, 511], [2, 3, 5, 9, 17, 32, 33, 48, 56, 60, 
    62, 63, 65, 80, 88, 92, 94, 95, 104, 108, 110, 111, 116, 118, 119, 122, 123, 125, 129, 144, 152, 156, 158, 159, 168, 172, 174, 175, 180, 182, 183, 186, 187, 189, 200, 204, 206, 207, 212, 214, 215, 218, 219, 221, 228, 230, 231, 234, 235, 237, 242, 243, 245, 249, 257, 272, 280, 284, 286, 287, 296, 300, 302, 303, 308, 310, 311, 314, 315, 317, 328, 332, 334, 335, 340, 342, 343, 346, 347, 349, 356, 358, 359, 362, 363, 365, 370, 371, 373, 377, 392, 396, 398, 399, 
    404, 406, 407, 410, 411, 413, 420, 422, 423, 426, 427, 429, 434, 435, 437, 441, 452, 454, 455, 458, 459, 461, 466, 467, 469, 473, 482, 483, 485, 489, 497, 512], [4, 6, 7, 10, 11, 13, 18, 19, 21, 25, 34, 35, 37, 41, 49, 64, 66, 67, 
    69, 73, 81, 96, 97, 112, 120, 124, 126, 127, 130, 131, 133, 137, 145, 160, 161, 176, 184, 188, 190, 191, 193, 208, 216, 220, 222, 223, 232, 236, 238, 239, 244, 246, 247, 250, 251, 253, 258, 259, 261, 265, 273, 288, 289, 304, 312, 316, 318, 319, 321, 336, 344, 348, 350, 351, 360, 364, 366, 367, 372, 374, 375, 378, 379, 381, 385, 400, 408, 412, 414, 415, 424, 428, 430, 431, 436, 438, 439, 442, 443, 445, 456, 460, 462, 463, 468, 470, 471, 474, 475, 477, 484, 486, 487, 490, 491, 493, 498, 499, 501, 505], [8, 12, 14, 15, 20, 22, 23, 26, 27, 29, 36, 38, 39, 42, 43, 45, 50, 51, 53, 57, 68, 70, 71, 74, 75, 77, 82, 83, 85, 89, 98, 99, 101, 105, 113, 128, 132, 134, 135, 138, 139, 141, 146, 147, 149, 153, 162, 163, 165, 169, 177, 192, 194, 195, 197, 201, 209, 224, 225, 240, 248, 252, 254, 255, 260, 262, 263, 266, 267, 269, 274, 275, 277, 281, 290, 291, 293, 297, 305, 320, 322, 323, 325, 329, 337, 352, 353, 368, 376, 380, 
    382, 383, 386, 387, 389, 393, 401, 416, 417, 432, 440, 444, 446, 447, 449, 464, 472, 476, 478, 479, 488, 492, 494, 495, 500, 502, 503, 506, 507, 509]]
    pos = [@variable(model, [1:length(blocks[i]), 1:length(blocks[i])], PSD) for i =1:4]
    for i = 0:3, j = 0:3, k = 0:3, l = 0:3, s = 0:3, t = 0:3, u = 0:3, v = 0:3, w = 0:3
        ind = [i,j,k,l,s,t,u,v,w]
        if all(x->iseven(sum(ind .== x)), 1:3)
            inx = ind .!= 0
            Locb = bfind(tsupp, ltsupp, reduce4(UInt16.(3*(Vector(1:9)[inx] .- 1) + ind[inx]), L, lattice=lattice))
            tp = real(kron(Pauli[i+1], Pauli[j+1], Pauli[k+1], Pauli[l+1], Pauli[s+1], Pauli[t+1], Pauli[u+1], Pauli[v+1], Pauli[w+1]))
            for p = 1:4
                @inbounds add_to_expression!(cons[Locb], sum(tp[blocks[p], blocks[p]].*pos[p]))
            end
        end
    end
end

function posepsd10!(model, cons, tsupp, L; lattice="chain")
    ltsupp = length(tsupp)
    Pauli = Matrix{Complex{Int8}}[[1 0; 0 1], [0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]
    blocks = [[1, 16, 24, 28, 30, 31, 40, 44, 46, 47, 52, 54, 55, 58, 59, 61, 72, 76, 78, 79, 84, 86, 87, 90, 91, 93, 100, 102, 103, 106, 107, 109, 114, 115, 117, 121, 136, 140, 142, 143, 148, 150, 151, 154, 155, 157, 164, 166, 167, 170, 171, 173, 178, 179, 181, 185, 196, 198, 199, 202, 203, 205, 210, 211, 213, 217, 226, 227, 229, 233, 241, 256, 264, 268, 270, 271, 276, 278, 279, 282, 283, 285, 292, 294, 295, 298, 299, 301, 306, 307, 309, 313, 324, 326, 327, 330, 331, 333, 338, 339, 341, 345, 354, 355, 357, 361, 369, 384, 388, 390, 391, 394, 395, 397, 402, 403, 405, 409, 418, 419, 421, 425, 433, 448, 450, 451, 453, 457, 465, 480, 481, 496, 504, 508, 510, 511, 520, 524, 526, 527, 532, 534, 535, 538, 539, 541, 548, 550, 551, 554, 555, 557, 562, 563, 565, 569, 580, 582, 583, 586, 587, 589, 594, 595, 597, 601, 610, 611, 613, 617, 625, 640, 644, 646, 647, 650, 651, 653, 658, 659, 661, 665, 674, 675, 677, 681, 689, 704, 706, 707, 709, 713, 721, 736, 737, 752, 760, 764, 766, 767, 772, 774, 775, 778, 779, 781, 786, 787, 789, 793, 802, 803, 805, 809, 817, 832, 834, 835, 837, 841, 849, 864, 865, 880, 888, 892, 894, 895, 898, 899, 901, 905, 913, 928, 929, 944, 
    952, 956, 958, 959, 961, 976, 984, 988, 990, 991, 1000, 1004, 1006, 1007, 1012, 1014, 1015, 1018, 1019, 1021], [2, 3, 5, 9, 17, 32, 33, 48, 56, 60, 62, 63, 65, 80, 88, 92, 94, 95, 104, 108, 110, 111, 116, 118, 119, 122, 123, 125, 129, 144, 152, 156, 158, 159, 168, 172, 174, 175, 180, 182, 183, 186, 187, 189, 200, 204, 206, 207, 212, 214, 215, 218, 219, 221, 228, 230, 231, 234, 235, 237, 242, 243, 245, 249, 257, 272, 280, 284, 286, 287, 296, 300, 302, 303, 308, 310, 311, 314, 315, 317, 328, 332, 334, 335, 340, 342, 343, 346, 347, 349, 356, 358, 359, 362, 363, 365, 370, 371, 373, 377, 392, 396, 398, 399, 404, 406, 407, 410, 411, 413, 420, 422, 423, 426, 427, 429, 434, 435, 437, 441, 452, 454, 455, 458, 459, 461, 466, 467, 469, 473, 482, 483, 485, 489, 497, 512, 513, 528, 536, 540, 542, 543, 552, 556, 558, 559, 564, 566, 567, 570, 571, 573, 584, 588, 590, 591, 596, 598, 599, 602, 603, 605, 612, 614, 615, 618, 619, 621, 626, 627, 629, 633, 648, 652, 654, 655, 660, 662, 663, 666, 667, 669, 676, 678, 679, 682, 683, 685, 690, 691, 693, 697, 708, 710, 711, 714, 715, 717, 722, 723, 725, 729, 738, 739, 741, 745, 753, 768, 776, 780, 782, 783, 788, 
    790, 791, 794, 795, 797, 804, 806, 807, 810, 811, 813, 818, 819, 821, 825, 836, 838, 839, 842, 843, 845, 850, 851, 853, 857, 866, 867, 869, 873, 881, 896, 900, 902, 903, 906, 907, 909, 914, 915, 917, 921, 930, 931, 933, 937, 945, 960, 962, 963, 965, 969, 977, 992, 993, 1008, 1016, 1020, 1022, 1023], [4, 6, 7, 10, 11, 13, 18, 19, 21, 25, 34, 35, 37, 41, 49, 64, 66, 67, 69, 73, 81, 96, 97, 112, 120, 124, 126, 127, 130, 131, 133, 137, 145, 160, 161, 176, 184, 188, 190, 191, 193, 208, 216, 220, 222, 223, 232, 236, 238, 239, 244, 246, 247, 250, 251, 253, 258, 259, 261, 265, 273, 288, 289, 304, 312, 316, 318, 319, 321, 336, 344, 348, 350, 351, 360, 364, 366, 367, 372, 374, 375, 378, 379, 381, 385, 400, 408, 412, 414, 415, 424, 428, 430, 431, 436, 438, 439, 442, 443, 445, 456, 460, 462, 463, 468, 470, 471, 474, 475, 477, 484, 486, 487, 490, 491, 493, 498, 499, 501, 505, 514, 515, 517, 521, 529, 544, 545, 560, 568, 572, 574, 575, 577, 592, 600, 604, 606, 607, 616, 620, 622, 623, 628, 630, 631, 634, 635, 637, 641, 656, 664, 668, 670, 671, 680, 684, 686, 687, 692, 694, 695, 698, 699, 701, 712, 716, 718, 719, 724, 726, 727, 730, 731, 733, 740, 742, 743, 746, 747, 749, 754, 755, 757, 761, 769, 784, 792, 796, 798, 799, 808, 812, 814, 815, 820, 822, 823, 826, 827, 829, 840, 844, 846, 847, 852, 854, 855, 858, 859, 861, 868, 870, 871, 874, 875, 877, 882, 883, 885, 889, 904, 908, 
    910, 911, 916, 918, 919, 922, 923, 925, 932, 934, 935, 938, 939, 941, 946, 947, 949, 953, 964, 966, 967, 970, 971, 973, 978, 979, 981, 985, 994, 995, 997, 1001, 1009, 1024], [8, 12, 14, 15, 20, 22, 23, 26, 27, 29, 36, 38, 39, 42, 43, 45, 50, 51, 53, 57, 68, 70, 71, 74, 75, 77, 82, 83, 85, 89, 98, 99, 101, 105, 113, 128, 132, 134, 135, 138, 139, 141, 146, 147, 149, 153, 162, 163, 165, 169, 177, 192, 194, 195, 197, 201, 209, 224, 225, 240, 248, 252, 254, 255, 
    260, 262, 263, 266, 267, 269, 274, 275, 277, 281, 290, 291, 293, 297, 305, 320, 322, 323, 325, 329, 337, 352, 353, 368, 376, 380, 382, 383, 386, 387, 389, 393, 401, 416, 417, 432, 440, 444, 446, 447, 449, 464, 472, 476, 478, 479, 488, 492, 494, 495, 500, 502, 503, 506, 507, 509, 516, 518, 519, 522, 523, 525, 530, 531, 533, 537, 546, 547, 549, 553, 561, 576, 578, 579, 581, 585, 593, 608, 609, 624, 632, 636, 638, 639, 642, 643, 645, 649, 657, 672, 673, 688, 696, 700, 702, 703, 705, 720, 728, 732, 734, 735, 744, 748, 750, 751, 756, 758, 759, 762, 763, 765, 770, 771, 773, 777, 785, 800, 801, 816, 824, 828, 830, 831, 833, 848, 856, 860, 862, 863, 872, 876, 878, 879, 884, 886, 887, 890, 891, 893, 897, 912, 920, 924, 926, 927, 936, 940, 942, 943, 948, 950, 951, 954, 955, 957, 968, 972, 974, 975, 980, 982, 983, 986, 987, 989, 996, 998, 999, 1002, 1003, 1005, 1010, 1011, 1013, 1017]]
    pos = [@variable(model, [1:length(blocks[i]), 1:length(blocks[i])], PSD) for i =1:4]
    for i = 0:3, j = 0:3, k = 0:3, l = 0:3, s = 0:3, t = 0:3, u = 0:3, v = 0:3, w = 0:3, z = 0:3
        ind = [i,j,k,l,s,t,u,v,w,z]
        if all(x->iseven(sum(ind .== x)), 1:3)
            inx = ind .!= 0
            Locb = bfind(tsupp, ltsupp, reduce4(UInt16.(3*(Vector(1:10)[inx] .- 1) + ind[inx]), L, lattice=lattice))
            tp = real(kron(Pauli[i+1], Pauli[j+1], Pauli[k+1], Pauli[l+1], Pauli[s+1], Pauli[t+1], Pauli[u+1], Pauli[v+1], Pauli[w+1], Pauli[z+1]))
            for p = 1:4
                @inbounds add_to_expression!(cons[Locb], sum(tp[blocks[p], blocks[p]].*pos[p]))
            end
        end
    end
end
