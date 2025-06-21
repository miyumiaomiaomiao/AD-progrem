library(Seurat)
library(harmony)


sample_1 <- Read10X("C:/Seurat_fold(MTG)/H18.30.001_L2_3_IT")
seurat_1 <- CreateSeuratObject(counts = sample_1, project = "H18.30.001_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_1), size = 0.2 * length(Cells(seurat_1)))
seurat_1_downsampled <- subset(seurat_1, cells = cells_to_keep)
rm(sample_1)
rm(seurat_1)


sample_2 <- Read10X("C:/Seurat_fold(MTG)/H18.30.002_L2_3_IT")
seurat_2 <- CreateSeuratObject(counts = sample_2, project = "H18.30.002_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_2), size = 0.2 * length(Cells(seurat_2)))
seurat_2_downsampled <- subset(seurat_2, cells = cells_to_keep)
rm(sample_2)
rm(seurat_2)


sample_3 <- Read10X("C:/Seurat_fold(MTG)/H19.30.001_L2_3_IT")
seurat_3 <- CreateSeuratObject(counts = sample_3, project = "H19.30.001_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_3), size = 0.2 * length(Cells(seurat_3)))
seurat_3_downsampled <- subset(seurat_3, cells = cells_to_keep)
rm(sample_3)
rm(seurat_3)


sample_4 <- Read10X("C:/Seurat_fold(MTG)/H19.30.002_L2_3_IT")
seurat_4 <- CreateSeuratObject(counts = sample_4, project = "H19.30.002_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_4), size = 0.2 * length(Cells(seurat_4)))
seurat_4_downsampled <- subset(seurat_4, cells = cells_to_keep)
rm(sample_4)
rm(seurat_4)


sample_5 <- Read10X("C:/Seurat_fold(MTG)/H19.33.004_L2_3_IT")
seurat_5 <- CreateSeuratObject(counts = sample_5, project = "H19.33.004_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_5), size = 0.2 * length(Cells(seurat_5)))
seurat_5_downsampled <- subset(seurat_5, cells = cells_to_keep)
rm(sample_5)
rm(seurat_5)


sample_6 <- Read10X("C:/Seurat_fold(MTG)/H20.33.001_L2_3_IT")
seurat_6 <- CreateSeuratObject(counts = sample_6, project = "H20.33.001_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_6), size = 0.2 * length(Cells(seurat_6)))
seurat_6_downsampled <- subset(seurat_6, cells = cells_to_keep)
rm(sample_6)
rm(seurat_6)


sample_7 <- Read10X("C:/Seurat_fold(MTG)/H20.33.002_L2_3_IT")
seurat_7 <- CreateSeuratObject(counts = sample_7, project = "H20.33.002_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_7), size = 0.2 * length(Cells(seurat_7)))
seurat_7_downsampled <- subset(seurat_7, cells = cells_to_keep)
rm(sample_7)
rm(seurat_7)


sample_8 <- Read10X("C:/Seurat_fold(MTG)/H20.33.004_L2_3_IT")
seurat_8 <- CreateSeuratObject(counts = sample_8, project = "H20.33.004_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_8), size = 0.2 * length(Cells(seurat_8)))
seurat_8_downsampled <- subset(seurat_8, cells = cells_to_keep)
rm(sample_8)
rm(seurat_8)


sample_9 <- Read10X("C:/Seurat_fold(MTG)/H20.33.005_L2_3_IT")
seurat_9 <- CreateSeuratObject(counts = sample_9, project = "H20.33.005_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_9), size = 0.2 * length(Cells(seurat_9)))
seurat_9_downsampled <- subset(seurat_9, cells = cells_to_keep)
rm(sample_9)
rm(seurat_9)


sample_10 <- Read10X("C:/Seurat_fold(MTG)/H20.33.008_L2_3_IT")
seurat_10 <- CreateSeuratObject(counts = sample_10, project = "H20.33.008_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_10), size = 0.2 * length(Cells(seurat_10)))
seurat_10_downsampled <- subset(seurat_10, cells = cells_to_keep)
rm(sample_10)
rm(seurat_10)


sample_11 <- Read10X("C:/Seurat_fold(MTG)/H20.33.011_L2_3_IT")
seurat_11 <- CreateSeuratObject(counts = sample_11, project = "H20.33.011_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_11), size = 0.2 * length(Cells(seurat_11)))
seurat_11_downsampled <- subset(seurat_11, cells = cells_to_keep)
rm(sample_11)
rm(seurat_11)


sample_12 <- Read10X("C:/Seurat_fold(MTG)/H20.33.012_L2_3_IT")
seurat_12 <- CreateSeuratObject(counts = sample_12, project = "H20.33.012_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_12), size = 0.2 * length(Cells(seurat_12)))
seurat_12_downsampled <- subset(seurat_12, cells = cells_to_keep)
rm(sample_12)
rm(seurat_12)


sample_13 <- Read10X("C:/Seurat_fold(MTG)/H20.33.013_L2_3_IT")
seurat_13 <- CreateSeuratObject(counts = sample_13, project = "H20.33.013_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_13), size = 0.2 * length(Cells(seurat_13)))
seurat_13_downsampled <- subset(seurat_13, cells = cells_to_keep)
rm(sample_13)
rm(seurat_13)


sample_14 <- Read10X("C:/Seurat_fold(MTG)/H20.33.014_L2_3_IT")
seurat_14 <- CreateSeuratObject(counts = sample_14, project = "H20.33.014_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_14), size = 0.2 * length(Cells(seurat_14)))
seurat_14_downsampled <- subset(seurat_14, cells = cells_to_keep)
rm(sample_14)
rm(seurat_14)


sample_15 <- Read10X("C:/Seurat_fold(MTG)/H20.33.015_L2_3_IT")
seurat_15 <- CreateSeuratObject(counts = sample_15, project = "H20.33.015_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_15), size = 0.2 * length(Cells(seurat_15)))
seurat_15_downsampled <- subset(seurat_15, cells = cells_to_keep)
rm(sample_15)
rm(seurat_15)


sample_16 <- Read10X("C:/Seurat_fold(MTG)/H20.33.016_L2_3_IT")
seurat_16 <- CreateSeuratObject(counts = sample_16, project = "H20.33.016_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_16), size = 0.2 * length(Cells(seurat_16)))
seurat_16_downsampled <- subset(seurat_16, cells = cells_to_keep)
rm(sample_16)
rm(seurat_16)

sample_17 <- Read10X("C:/Seurat_fold(MTG)/H20.33.017_L2_3_IT")
seurat_17 <- CreateSeuratObject(counts = sample_17, project = "H20.33.017_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_17), size = 0.2 * length(Cells(seurat_17)))
seurat_17_downsampled <- subset(seurat_17, cells = cells_to_keep)
rm(sample_17)
rm(seurat_17)


sample_18 <- Read10X("C:/Seurat_fold(MTG)/H20.33.018_L2_3_IT")
seurat_18 <- CreateSeuratObject(counts = sample_18, project = "H20.33.018_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_18), size = 0.2 * length(Cells(seurat_18)))
seurat_18_downsampled <- subset(seurat_18, cells = cells_to_keep)
rm(sample_18)
rm(seurat_18)


sample_19 <- Read10X("C:/Seurat_fold(MTG)/H20.33.019_L2_3_IT")
seurat_19 <- CreateSeuratObject(counts = sample_19, project = "H20.33.019_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_19), size = 0.2 * length(Cells(seurat_19)))
seurat_19_downsampled <- subset(seurat_19, cells = cells_to_keep)
rm(sample_19)
rm(seurat_19)


sample_20 <- Read10X("C:/Seurat_fold(MTG)/H20.33.020_L2_3_IT")
seurat_20 <- CreateSeuratObject(counts = sample_20, project = "H20.33.020_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_20), size = 0.2 * length(Cells(seurat_20)))
seurat_20_downsampled <- subset(seurat_20, cells = cells_to_keep)
rm(sample_20)
rm(seurat_20)


sample_21 <- Read10X("C:/Seurat_fold(MTG)/H20.33.024_L2_3_IT")
seurat_21 <- CreateSeuratObject(counts = sample_21, project = "H20.33.024_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_21), size = 0.2 * length(Cells(seurat_21)))
seurat_21_downsampled <- subset(seurat_21, cells = cells_to_keep)
rm(sample_21)
rm(seurat_21)

sample_22 <- Read10X("C:/Seurat_fold(MTG)/H20.33.025_L2_3_IT")
seurat_22 <- CreateSeuratObject(counts = sample_22, project = "H20.33.025_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_22), size = 0.2 * length(Cells(seurat_22)))
seurat_22_downsampled <- subset(seurat_22, cells = cells_to_keep)
rm(sample_22)
rm(seurat_22)


sample_23 <- Read10X("C:/Seurat_fold(MTG)/H20.33.026_L2_3_IT")
seurat_23 <- CreateSeuratObject(counts = sample_23, project = "H20.33.026_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_23), size = 0.2 * length(Cells(seurat_23)))
seurat_23_downsampled <- subset(seurat_23, cells = cells_to_keep)
rm(sample_23)
rm(seurat_23)


sample_24 <- Read10X("C:/Seurat_fold(MTG)/H20.33.027_L2_3_IT")
seurat_24 <- CreateSeuratObject(counts = sample_24, project = "H20.33.027_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_24), size = 0.2 * length(Cells(seurat_24)))
seurat_24_downsampled <- subset(seurat_24, cells = cells_to_keep)
rm(sample_24)
rm(seurat_24)


sample_26 <- Read10X("C:/Seurat_fold(MTG)/H20.33.028_L2_3_IT")
seurat_26 <- CreateSeuratObject(counts = sample_26, project = "H20.33.028_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_26), size = 0.2 * length(Cells(seurat_26)))
seurat_26_downsampled <- subset(seurat_26, cells = cells_to_keep)
rm(sample_26)
rm(seurat_26)


sample_27 <- Read10X("C:/Seurat_fold(MTG)/H20.33.029_L2_3_IT")
seurat_27 <- CreateSeuratObject(counts = sample_27, project = "H20.33.029_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_27), size = 0.2 * length(Cells(seurat_27)))
seurat_27_downsampled <- subset(seurat_27, cells = cells_to_keep)
rm(sample_27)
rm(seurat_27)


sample_29 <- Read10X("C:/Seurat_fold(MTG)/H20.33.030_L2_3_IT")
seurat_29 <- CreateSeuratObject(counts = sample_29, project = "H20.33.030_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_29), size = 0.2 * length(Cells(seurat_29)))
seurat_29_downsampled <- subset(seurat_29, cells = cells_to_keep)
rm(sample_29)
rm(seurat_29)


sample_30 <- Read10X("C:/Seurat_fold(MTG)/H20.33.031_L2_3_IT")
seurat_30 <- CreateSeuratObject(counts = sample_30, project = "H20.33.031_L2_3_IT")  # 修正 project 后缀
cells_to_keep <- sample(Cells(seurat_30), size = 0.2 * length(Cells(seurat_30)))
seurat_30_downsampled <- subset(seurat_30, cells = cells_to_keep)
rm(sample_30)


sample_31 <- Read10X("C:/Seurat_fold(MTG)/H20.33.032_L2_3_IT")
seurat_31 <- CreateSeuratObject(counts = sample_31, project = "H20.33.032_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_31), size = 0.2 * length(Cells(seurat_31)))
seurat_31_downsampled <- subset(seurat_31, cells = cells_to_keep)
rm(sample_31)
rm(seurat_31)

sample_32 <- Read10X("C:/Seurat_fold(MTG)/H20.33.033_L2_3_IT")
seurat_32 <- CreateSeuratObject(counts = sample_32, project = "H20.33.033_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_32), size = 0.2 * length(Cells(seurat_32)))
seurat_32_downsampled <- subset(seurat_32, cells = cells_to_keep)
rm(sample_32)
rm(seurat_32)


sample_33 <- Read10X("C:/Seurat_fold(MTG)/H20.33.034_L2_3_IT")
seurat_33 <- CreateSeuratObject(counts = sample_33, project = "H20.33.034_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_33), size = 0.2 * length(Cells(seurat_33)))
seurat_33_downsampled <- subset(seurat_33, cells = cells_to_keep)
rm(sample_33)
rm(seurat_33)


sample_34 <- Read10X("C:/Seurat_fold(MTG)/H20.33.035_L2_3_IT")
seurat_34 <- CreateSeuratObject(counts = sample_34, project = "H20.33.035_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_34), size = 0.2 * length(Cells(seurat_34)))
seurat_34_downsampled <- subset(seurat_34, cells = cells_to_keep)
rm(sample_34)
rm(seurat_34)


sample_35 <- Read10X("C:/Seurat_fold(MTG)/H20.33.036_L2_3_IT")
seurat_35 <- CreateSeuratObject(counts = sample_35, project = "H20.33.036_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_35), size = 0.2 * length(Cells(seurat_35)))
seurat_35_downsampled <- subset(seurat_35, cells = cells_to_keep)
rm(sample_35)
rm(seurat_35)


sample_36 <- Read10X("C:/Seurat_fold(MTG)/H20.33.037_L2_3_IT")
seurat_36 <- CreateSeuratObject(counts = sample_36, project = "H20.33.037_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_36), size = 0.2 * length(Cells(seurat_36)))
seurat_36_downsampled <- subset(seurat_36, cells = cells_to_keep)
rm(sample_36)
rm(seurat_36)


sample_37 <- Read10X("C:/Seurat_fold(MTG)/H20.33.038_L2_3_IT")
seurat_37 <- CreateSeuratObject(counts = sample_37, project = "H20.33.038_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_37), size = 0.2 * length(Cells(seurat_37)))
seurat_37_downsampled <- subset(seurat_37, cells = cells_to_keep)
rm(sample_37)
rm(seurat_37)


sample_38 <- Read10X("C:/Seurat_fold(MTG)/H20.33.039_L2_3_IT")
seurat_38 <- CreateSeuratObject(counts = sample_38, project = "H20.33.039_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_38), size = 0.2 * length(Cells(seurat_38)))
seurat_38_downsampled <- subset(seurat_38, cells = cells_to_keep)
rm(sample_38)
rm(seurat_38)

sample_39 <- Read10X("C:/Seurat_fold(MTG)/H20.33.040_L2_3_IT")
seurat_39 <- CreateSeuratObject(counts = sample_39, project = "H20.33.040_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_39), size = 0.2 * length(Cells(seurat_39)))
seurat_39_downsampled <- subset(seurat_39, cells = cells_to_keep)
rm(sample_39)
rm(seurat_39)


sample_40 <- Read10X("C:/Seurat_fold(MTG)/H20.33.041_L2_3_IT")
seurat_40 <- CreateSeuratObject(counts = sample_40, project = "H20.33.041_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_40), size = 0.2 * length(Cells(seurat_40)))
seurat_40_downsampled <- subset(seurat_40, cells = cells_to_keep)
rm(sample_40)
rm(seurat_40)


sample_41 <- Read10X("C:/Seurat_fold(MTG)/H20.33.043_L2_3_IT")
seurat_41 <- CreateSeuratObject(counts = sample_41, project = "H20.33.043_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_41), size = 0.2 * length(Cells(seurat_41)))
seurat_41_downsampled <- subset(seurat_41, cells = cells_to_keep)
rm(sample_41)
rm(seurat_41)


sample_42 <- Read10X("C:/Seurat_fold(MTG)/H20.33.044_L2_3_IT")
seurat_42 <- CreateSeuratObject(counts = sample_42, project = "H20.33.044_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_42), size = 0.2 * length(Cells(seurat_42)))
seurat_42_downsampled <- subset(seurat_42, cells = cells_to_keep)
rm(sample_42)
rm(seurat_42)

sample_43 <- Read10X("C:/Seurat_fold(MTG)/H20.33.045_L2_3_IT")
seurat_43 <- CreateSeuratObject(counts = sample_43, project = "H20.33.045_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_43), size = 0.2 * length(Cells(seurat_43)))
seurat_43_downsampled <- subset(seurat_43, cells = cells_to_keep)
rm(sample_43)
rm(seurat_43)


sample_44 <- Read10X("C:/Seurat_fold(MTG)/H20.33.046_L2_3_IT")
seurat_44 <- CreateSeuratObject(counts = sample_44, project = "H20.33.046_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_44), size = 0.2 * length(Cells(seurat_44)))
seurat_44_downsampled <- subset(seurat_44, cells = cells_to_keep)
rm(sample_44)
rm(seurat_44)

sample_45 <- Read10X("C:/Seurat_fold(MTG)/H21.33.001_L2_3_IT")
seurat_45 <- CreateSeuratObject(counts = sample_45, project = "H21.33.001_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_45), size = 0.2 * length(Cells(seurat_45)))
seurat_45_downsampled <- subset(seurat_45, cells = cells_to_keep)
rm(sample_45)
rm(seurat_45)


sample_46 <- Read10X("C:/Seurat_fold(MTG)/H21.33.002_L2_3_IT")
seurat_46 <- CreateSeuratObject(counts = sample_46, project = "H21.33.002_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_46), size = 0.2 * length(Cells(seurat_46)))
seurat_46_downsampled <- subset(seurat_46, cells = cells_to_keep)
rm(sample_46)
rm(seurat_46)


sample_47 <- Read10X("C:/Seurat_fold(MTG)/H21.33.003_L2_3_IT")
seurat_47 <- CreateSeuratObject(counts = sample_47, project = "H21.33.003_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_47), size = 0.2 * length(Cells(seurat_47)))
seurat_47_downsampled <- subset(seurat_47, cells = cells_to_keep)
rm(sample_47)
rm(seurat_47)


sample_48 <- Read10X("C:/Seurat_fold(MTG)/H21.33.004_L2_3_IT")
seurat_48 <- CreateSeuratObject(counts = sample_48, project = "H21.33.004_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_48), size = 0.2 * length(Cells(seurat_48)))
seurat_48_downsampled <- subset(seurat_48, cells = cells_to_keep)
rm(sample_48)
rm(seurat_48)


sample_49 <- Read10X("C:/Seurat_fold(MTG)/H21.33.005_L2_3_IT")
seurat_49 <- CreateSeuratObject(counts = sample_49, project = "H21.33.005_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_49), size = 0.2 * length(Cells(seurat_49)))
seurat_49_downsampled <- subset(seurat_49, cells = cells_to_keep)
rm(sample_49)
rm(seurat_49)


sample_50 <- Read10X("C:/Seurat_fold(MTG)/H21.33.006_L2_3_IT")
seurat_50 <- CreateSeuratObject(counts = sample_50, project = "H21.33.006_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_50), size = 0.2 * length(Cells(seurat_50)))
seurat_50_downsampled <- subset(seurat_50, cells = cells_to_keep)
rm(sample_50)
rm(seurat_50)


sample_51 <- Read10X("C:/Seurat_fold(MTG)/H21.33.007_L2_3_IT")
seurat_51 <- CreateSeuratObject(counts = sample_51, project = "H21.33.007_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_51), size = 0.2 * length(Cells(seurat_51)))
seurat_51_downsampled <- subset(seurat_51, cells = cells_to_keep)
rm(sample_51)
rm(seurat_51)


sample_52 <- Read10X("C:/Seurat_fold(MTG)/H21.33.008_L2_3_IT")
seurat_52 <- CreateSeuratObject(counts = sample_52, project = "H21.33.008_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_52), size = 0.2 * length(Cells(seurat_52)))
seurat_52_downsampled <- subset(seurat_52, cells = cells_to_keep)
rm(sample_52)
rm(seurat_52)


sample_53 <- Read10X("C:/Seurat_fold(MTG)/H21.33.009_L2_3_IT")
seurat_53 <- CreateSeuratObject(counts = sample_53, project = "H21.33.009_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_53), size = 0.2 * length(Cells(seurat_53)))
seurat_53_downsampled <- subset(seurat_53, cells = cells_to_keep)
rm(sample_53)
rm(seurat_53)
sample_54 <- Read10X("C:/Seurat_fold(MTG)/H21.33.010_L2_3_IT")
seurat_54 <- CreateSeuratObject(counts = sample_54, project = "H21.33.010_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_54), size = 0.2 * length(Cells(seurat_54)))
seurat_54_downsampled <- subset(seurat_54, cells = cells_to_keep)
rm(sample_54)
rm(seurat_54)


sample_55 <- Read10X("C:/Seurat_fold(MTG)/H21.33.011_L2_3_IT")
seurat_55 <- CreateSeuratObject(counts = sample_55, project = "H21.33.011_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_55), size = 0.2 * length(Cells(seurat_55)))
seurat_55_downsampled <- subset(seurat_55, cells = cells_to_keep)
rm(sample_55)
rm(seurat_55)


sample_56 <- Read10X("C:/Seurat_fold(MTG)/H21.33.012_L2_3_IT")
seurat_56 <- CreateSeuratObject(counts = sample_56, project = "H21.33.012_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_56), size = 0.2 * length(Cells(seurat_56)))
seurat_56_downsampled <- subset(seurat_56, cells = cells_to_keep)
rm(sample_56)
rm(seurat_56)

sample_57 <- Read10X("C:/Seurat_fold(MTG)/H21.33.013_L2_3_IT")
seurat_57 <- CreateSeuratObject(counts = sample_57, project = "H21.33.013_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_57), size = 0.2 * length(Cells(seurat_57)))
seurat_57_downsampled <- subset(seurat_57, cells = cells_to_keep)
rm(sample_57)
rm(seurat_57)


sample_58 <- Read10X("C:/Seurat_fold(MTG)/H21.33.014_L2_3_IT")
seurat_58 <- CreateSeuratObject(counts = sample_58, project = "H21.33.014_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_58), size = 0.2 * length(Cells(seurat_58)))
seurat_58_downsampled <- subset(seurat_58, cells = cells_to_keep)
rm(sample_58)
rm(seurat_58)

sample_59 <- Read10X("C:/Seurat_fold(MTG)/H21.33.015_L2_3_IT")
seurat_59 <- CreateSeuratObject(counts = sample_59, project = "H21.33.015_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_59), size = 0.2 * length(Cells(seurat_59)))
seurat_59_downsampled <- subset(seurat_59, cells = cells_to_keep)
rm(sample_59)
rm(seurat_59)

sample_60 <- Read10X("C:/Seurat_fold(MTG)/H21.33.016_L2_3_IT")
seurat_60 <- CreateSeuratObject(counts = sample_60, project = "H21.33.016_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_60), size = 0.2 * length(Cells(seurat_60)))
seurat_60_downsampled <- subset(seurat_60, cells = cells_to_keep)
rm(sample_60)
rm(seurat_60)

sample_61 <- Read10X("C:/Seurat_fold(MTG)/H21.33.017_L2_3_IT")
seurat_61 <- CreateSeuratObject(counts = sample_61, project = "H21.33.017_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_61), size = 0.2 * length(Cells(seurat_61)))
seurat_61_downsampled <- subset(seurat_61, cells = cells_to_keep)
rm(sample_61)
rm(seurat_61)

sample_62 <- Read10X("C:/Seurat_fold(MTG)/H21.33.018_L2_3_IT")
seurat_62 <- CreateSeuratObject(counts = sample_62, project = "H21.33.018_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_62), size = 0.2 * length(Cells(seurat_62)))
seurat_62_downsampled <- subset(seurat_62, cells = cells_to_keep)
rm(sample_62)

sample_63 <- Read10X("C:/Seurat_fold(MTG)/H21.33.019_L2_3_IT")
seurat_63 <- CreateSeuratObject(counts = sample_63, project = "H21.33.019_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_63), size = 0.2 * length(Cells(seurat_63)))
seurat_63_downsampled <- subset(seurat_63, cells = cells_to_keep)
rm(sample_63)
rm(seurat_63)

sample_64 <- Read10X("C:/Seurat_fold(MTG)/H21.33.020_L2_3_IT")
seurat_64 <- CreateSeuratObject(counts = sample_64, project = "H21.33.020_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_64), size = 0.2 * length(Cells(seurat_64)))
seurat_64_downsampled <- subset(seurat_64, cells = cells_to_keep)
rm(sample_64)
rm(seurat_64)


sample_65 <- Read10X("C:/Seurat_fold(MTG)/H21.33.021_L2_3_IT")
seurat_65 <- CreateSeuratObject(counts = sample_65, project = "H21.33.021_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_65), size = 0.2 * length(Cells(seurat_65)))
seurat_65_downsampled <- subset(seurat_65, cells = cells_to_keep)
rm(sample_65)
rm(seurat_65)

sample_66 <- Read10X("C:/Seurat_fold(MTG)/H21.33.022_L2_3_IT")
seurat_66 <- CreateSeuratObject(counts = sample_66, project = "H21.33.022_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_66), size = 0.2 * length(Cells(seurat_66)))
seurat_66_downsampled <- subset(seurat_66, cells = cells_to_keep)
rm(sample_66)
rm(seurat_66)


sample_67 <- Read10X("C:/Seurat_fold(MTG)/H21.33.023_L2_3_IT")
seurat_67 <- CreateSeuratObject(counts = sample_67, project = "H21.33.023_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_67), size = 0.2 * length(Cells(seurat_67)))
seurat_67_downsampled <- subset(seurat_67, cells = cells_to_keep)
rm(sample_67)
rm(seurat_67)

sample_68 <- Read10X("C:/Seurat_fold(MTG)/H21.33.025_L2_3_IT")
seurat_68 <- CreateSeuratObject(counts = sample_68, project = "H21.33.025_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_68), size = 0.2 * length(Cells(seurat_68)))
seurat_68_downsampled <- subset(seurat_68, cells = cells_to_keep)
rm(sample_68)
rm(seurat_68)


sample_69 <- Read10X("C:/Seurat_fold(MTG)/H21.33.026_L2_3_IT")
seurat_69 <- CreateSeuratObject(counts = sample_69, project = "H21.33.026_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_69), size = 0.2 * length(Cells(seurat_69)))
seurat_69_downsampled <- subset(seurat_69, cells = cells_to_keep)
rm(sample_69)
rm(seurat_69)


sample_70 <- Read10X("C:/Seurat_fold(MTG)/H21.33.027_L2_3_IT")
seurat_70 <- CreateSeuratObject(counts = sample_70, project = "H21.33.027_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_70), size = 0.2 * length(Cells(seurat_70)))
seurat_70_downsampled <- subset(seurat_70, cells = cells_to_keep)
rm(sample_70)
rm(seurat_70)


sample_71 <- Read10X("C:/Seurat_fold(MTG)/H21.33.028_L2_3_IT")
seurat_71 <- CreateSeuratObject(counts = sample_71, project = "H21.33.028_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_71), size = 0.2 * length(Cells(seurat_71)))
seurat_71_downsampled <- subset(seurat_71, cells = cells_to_keep)
rm(sample_71)
rm(seurat_71)


sample_72 <- Read10X("C:/Seurat_fold(MTG)/H21.33.029_L2_3_IT")
seurat_72 <- CreateSeuratObject(counts = sample_72, project = "H21.33.029_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_72), size = 0.2 * length(Cells(seurat_72)))
seurat_72_downsampled <- subset(seurat_72, cells = cells_to_keep)
rm(sample_72)
rm(seurat_72)


sample_73 <- Read10X("C:/Seurat_fold(MTG)/H21.33.030_L2_3_IT")
seurat_73 <- CreateSeuratObject(counts = sample_73, project = "H21.33.030_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_73), size = 0.2 * length(Cells(seurat_73)))
seurat_73_downsampled <- subset(seurat_73, cells = cells_to_keep)
rm(sample_73)
rm(seurat_73)


sample_74 <- Read10X("C:/Seurat_fold(MTG)/H21.33.031_L2_3_IT")
seurat_74 <- CreateSeuratObject(counts = sample_74, project = "H21.33.031_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_74), size = 0.2 * length(Cells(seurat_74)))
seurat_74_downsampled <- subset(seurat_74, cells = cells_to_keep)
rm(sample_74)
rm(seurat_74)


sample_75 <- Read10X("C:/Seurat_fold(MTG)/H21.33.032_L2_3_IT")
seurat_75 <- CreateSeuratObject(counts = sample_75, project = "H21.33.032_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_75), size = 0.2 * length(Cells(seurat_75)))
seurat_75_downsampled <- subset(seurat_75, cells = cells_to_keep)
rm(sample_75)
rm(seurat_75)


sample_76 <- Read10X("C:/Seurat_fold(MTG)/H21.33.033_L2_3_IT")
seurat_76 <- CreateSeuratObject(counts = sample_76, project = "H21.33.033_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_76), size = 0.2 * length(Cells(seurat_76)))
seurat_76_downsampled <- subset(seurat_76, cells = cells_to_keep)
rm(sample_76)
rm(seurat_76)


sample_77 <- Read10X("C:/Seurat_fold(MTG)/H21.33.034_L2_3_IT")
seurat_77 <- CreateSeuratObject(counts = sample_77, project = "H21.33.034_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_77), size = 0.2 * length(Cells(seurat_77)))
seurat_77_downsampled <- subset(seurat_77, cells = cells_to_keep)
rm(sample_77)
rm(seurat_77)


sample_78 <- Read10X("C:/Seurat_fold(MTG)/H21.33.035_L2_3_IT")
seurat_78 <- CreateSeuratObject(counts = sample_78, project = "H21.33.035_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_78), size = 0.2 * length(Cells(seurat_78)))
seurat_78_downsampled <- subset(seurat_78, cells = cells_to_keep)
rm(sample_78)
rm(seurat_78)

sample_79 <- Read10X("C:/Seurat_fold(MTG)/H21.33.036_L2_3_IT")
seurat_79 <- CreateSeuratObject(counts = sample_79, project = "H21.33.036_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_79), size = 0.2 * length(Cells(seurat_79)))
seurat_79_downsampled <- subset(seurat_79, cells = cells_to_keep)
rm(sample_79)
rm(seurat_79)


sample_81 <- Read10X("C:/Seurat_fold(MTG)/H21.33.037_L2_3_IT")
seurat_81 <- CreateSeuratObject(counts = sample_81, project = "H21.33.037_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_81), size = 0.2 * length(Cells(seurat_81)))
seurat_81_downsampled <- subset(seurat_81, cells = cells_to_keep)
rm(sample_81)
rm(seurat_81)


sample_82 <- Read10X("C:/Seurat_fold(MTG)/H21.33.038_L2_3_IT")
seurat_82 <- CreateSeuratObject(counts = sample_82, project = "H21.33.038_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_82), size = 0.2 * length(Cells(seurat_82)))
seurat_82_downsampled <- subset(seurat_82, cells = cells_to_keep)
rm(sample_82)
rm(seurat_82)

sample_83 <- Read10X("C:/Seurat_fold(MTG)/H21.33.039_L2_3_IT")
seurat_83 <- CreateSeuratObject(counts = sample_83, project = "H21.33.039_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_83), size = 0.2 * length(Cells(seurat_83)))
seurat_83_downsampled <- subset(seurat_83, cells = cells_to_keep)
rm(sample_83)
rm(seurat_83)

sample_84 <- Read10X("C:/Seurat_fold(MTG)/H21.33.040_L2_3_IT")
seurat_84 <- CreateSeuratObject(counts = sample_84, project = "H21.33.040_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_84), size = 0.2 * length(Cells(seurat_84)))
seurat_84_downsampled <- subset(seurat_84, cells = cells_to_keep)
rm(sample_84)
rm(seurat_84)

sample_85 <- Read10X("C:/Seurat_fold(MTG)/H21.33.041_L2_3_IT")
seurat_85 <- CreateSeuratObject(counts = sample_85, project = "H21.33.041_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_85), size = 0.2 * length(Cells(seurat_85)))
seurat_85_downsampled <- subset(seurat_85, cells = cells_to_keep)
rm(sample_85)
rm(seurat_85)

sample_86 <- Read10X("C:/Seurat_fold(MTG)/H21.33.042_L2_3_IT")
seurat_86 <- CreateSeuratObject(counts = sample_86, project = "H21.33.042_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_86), size = 0.2 * length(Cells(seurat_86)))
seurat_86_downsampled <- subset(seurat_86, cells = cells_to_keep)
rm(sample_86)
rm(seurat_86)

sample_87 <- Read10X("C:/Seurat_fold(MTG)/H21.33.043_L2_3_IT")
seurat_87 <- CreateSeuratObject(counts = sample_87, project = "H21.33.043_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_87), size = 0.2 * length(Cells(seurat_87)))
seurat_87_downsampled <- subset(seurat_87, cells = cells_to_keep)
rm(sample_87)
rm(seurat_87)

sample_88 <- Read10X("C:/Seurat_fold(MTG)/H21.33.044_L2_3_IT")
seurat_88 <- CreateSeuratObject(counts = sample_88, project = "H21.33.044_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_88), size = 0.2 * length(Cells(seurat_88)))
seurat_88_downsampled <- subset(seurat_88, cells = cells_to_keep)
rm(sample_88)
rm(seurat_88)

sample_89 <- Read10X("C:/Seurat_fold(MTG)/H21.33.045_L2_3_IT")
seurat_89 <- CreateSeuratObject(counts = sample_89, project = "H21.33.045_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_89), size = 0.2 * length(Cells(seurat_89)))
seurat_89_downsampled <- subset(seurat_89, cells = cells_to_keep)
rm(sample_89)
rm(seurat_89)

sample_90 <- Read10X("C:/Seurat_fold(MTG)/H21.33.046_L2_3_IT")
seurat_90 <- CreateSeuratObject(counts = sample_90, project = "H21.33.046_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_90), size = 0.2 * length(Cells(seurat_90)))
seurat_90_downsampled <- subset(seurat_90, cells = cells_to_keep)
rm(sample_90)
rm(seurat_90)

sample_91 <- Read10X("C:/Seurat_fold(MTG)/H21.33.047_L2_3_IT")
seurat_91 <- CreateSeuratObject(counts = sample_91, project = "H21.33.047_L2_3_IT")
cells_to_keep <- sample(Cells(seurat_91), size = 0.2 * length(Cells(seurat_91)))
seurat_91_downsampled <- subset(seurat_91, cells = cells_to_keep)
rm(sample_91)
rm(seurat_91)





dementia_merged <- merge(
  seurat_8_downsampled, 
  y = c(
    seurat_8_downsampled,seurat_11_downsampled,seurat_15_downsampled,
    seurat_16_downsampled,seurat_17_downsampled,seurat_18_downsampled,
    seurat_20_downsampled,seurat_23_downsampled,seurat_26_downsampled,
    seurat_27_downsampled,seurat_30_downsampled,seurat_32_downsampled,
    seurat_36_downsampled,seurat_37_downsampled,seurat_39_downsampled,
    seurat_40_downsampled,seurat_43_downsampled,seurat_44_downsampled,
    seurat_45_downsampled,seurat_46_downsampled,seurat_49_downsampled,
    seurat_51_downsampled,seurat_52_downsampled,seurat_53_downsampled,
    seurat_54_downsampled,seurat_56_downsampled,seurat_57_downsampled,
    seurat_60_downsampled,seurat_61_downsampled,seurat_62_downsampled,
    seurat_64_downsampled,seurat_65_downsampled,seurat_70_downsampled,
    seurat_72_downsampled,seurat_74_downsampled,seurat_77_downsampled,
    seurat_83_downsampled,seurat_86_downsampled,seurat_87_downsampled,
    seurat_88_downsampled,seurat_89_downsampled,seurat_90_downsampled
    
  ),
  add.cell.ids = paste0("Dementia_",1:43),
  project = "DementiaIntegrated")


no_dementia_merged <- merge(
  seurat_5_downsampled, 
  y = c(
    seurat_6_downsampled,seurat_7_downsampled,
    seurat_9_downsampled,seurat_10_downsampled,seurat_12_downsampled,
    seurat_13_downsampled,seurat_14_downsampled,seurat_19_downsampled,
    seurat_21_downsampled,seurat_22_downsampled,seurat_24_downsampled,
    seurat_29_downsampled,seurat_31_downsampled,seurat_33_downsampled,
    seurat_34_downsampled,seurat_35_downsampled,seurat_38_downsampled,
    seurat_41_downsampled,seurat_42_downsampled,seurat_47_downsampled,
    seurat_48_downsampled,seurat_50_downsampled,seurat_55_downsampled,
    seurat_58_downsampled,seurat_59_downsampled,seurat_63_downsampled,
    seurat_66_downsampled,seurat_67_downsampled,seurat_68_downsampled,
    seurat_69_downsampled,seurat_71_downsampled,seurat_73_downsampled,
    seurat_75_downsampled,seurat_76_downsampled,seurat_78_downsampled,
    seurat_79_downsampled,seurat_81_downsampled,seurat_82_downsampled,
    seurat_84_downsampled,seurat_85_downsampled,seurat_91_downsampled
  ),
  add.cell.ids = paste0("NoDementia_",1:42),
  project = "NoDementiaIntegrated")


combined_data <- merge(dementia_merged, y = no_dementia_merged, 
                       add.cell.ids = c("Dementia", "NoDementia"), 
                       project = "CombinedData")

combined_data <- NormalizeData(combined_data)
combined_data <- FindVariableFeatures(combined_data)
combined_data <- ScaleData(combined_data)
combined_data <- RunPCA(combined_data,reduction.name = "pca")
combined_data <- FindNeighbors(combined_data , dims = 1:20, reduction = "pca")
combined_data <- FindClusters(combined_data , resolution = 2, cluster.name = "dementiaornot_unintegrated_clusters")
combined_data <- RunUMAP(combined_data, dims = 1:20, reduction = "pca", reduction.name = "dementiaornot_unintegrated_umap")


DimPlot(combined_data, reduction = "dementiaornot_unintegrated_umap", group.by = "orig.ident")

class(combined_data)

combined_data <- IntegrateLayers(
  object = combined_data,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "combined_harmony",
  verbose = TRUE
)

combined_data <- FindNeighbors(combined_data, reduction = "combined_harmony", dims = 1:30)
combined_data <- FindClusters(combined_data, resolution = 2, cluster.name = "combined_harmony_clusters")
combined_data <- RunUMAP(combined_data, reduction = "combined_harmony", dims = 1:30, reduction.name = "combined_harmony_umap")


DimPlot(combined_data, reduction = "combined_harmony", group.by = "orig.ident")


conditions <- c(
  rep("no dementia", length(Cells(seurat_5_downsampled))),
  rep("no dementia", length(Cells(seurat_6_downsampled))),
  rep("no dementia", length(Cells(seurat_7_downsampled))),
  rep("dementia", length(Cells(seurat_8_downsampled))),
  rep("no dementia", length(Cells(seurat_9_downsampled))),
  rep("no dementia", length(Cells(seurat_10_downsampled))),
  rep("dementia", length(Cells(seurat_11_downsampled))),
  rep("no dementia", length(Cells(seurat_12_downsampled))),
  rep("no dementia", length(Cells(seurat_13_downsampled))),
  rep("no dementia", length(Cells(seurat_14_downsampled))),
  rep("dementia", length(Cells(seurat_15_downsampled))),
  rep("no dementia", length(Cells(seurat_16_downsampled))),
  rep("dementia", length(Cells(seurat_17_downsampled))),
  rep("no dementia", length(Cells(seurat_18_downsampled))),
  rep("no dementia", length(Cells(seurat_19_downsampled))),
  rep("dementia", length(Cells(seurat_20_downsampled))),
  rep("no dementia", length(Cells(seurat_21_downsampled))),
  rep("no dementia", length(Cells(seurat_22_downsampled))),
  rep("dementia", length(Cells(seurat_23_downsampled))),
  rep("no dementia", length(Cells(seurat_24_downsampled))),
  rep("dementia", length(Cells(seurat_26_downsampled))),
  rep("dementia", length(Cells(seurat_27_downsampled))),
  rep("no dementia", length(Cells(seurat_29_downsampled))),
  rep("dementia", length(Cells(seurat_30_downsampled))),
  rep("dementia", length(Cells(seurat_31_downsampled))),
  rep("dementia", length(Cells(seurat_32_downsampled))),
  rep("no dementia", length(Cells(seurat_33_downsampled))),
  rep("dementia", length(Cells(seurat_34_downsampled))),
  rep("dementia", length(Cells(seurat_35_downsampled))),
  rep("dementia", length(Cells(seurat_36_downsampled))),
  rep("dementia", length(Cells(seurat_37_downsampled))),
  rep("dementia", length(Cells(seurat_38_downsampled))),
  rep("dementia", length(Cells(seurat_39_downsampled))),
  rep("dementia", length(Cells(seurat_40_downsampled))),
  rep("dementia", length(Cells(seurat_41_downsampled))),
  rep("dementia", length(Cells(seurat_42_downsampled))),
  rep("dementia", length(Cells(seurat_43_downsampled))),
  rep("dementia", length(Cells(seurat_44_downsampled))),
  rep("dementia", length(Cells(seurat_45_downsampled))),
  rep("dementia", length(Cells(seurat_46_downsampled))),
  rep("dementia", length(Cells(seurat_47_downsampled))),
  rep("dementia", length(Cells(seurat_48_downsampled))),
  rep("no dementia", length(Cells(seurat_49_downsampled))),
  rep("dementia", length(Cells(seurat_50_downsampled))),
  rep("no dementia", length(Cells(seurat_51_downsampled))),
  rep("no dementia", length(Cells(seurat_52_downsampled))),
  rep("dementia", length(Cells(seurat_53_downsampled))),
  rep("no dementia", length(Cells(seurat_54_downsampled))),
  rep("dementia", length(Cells(seurat_55_downsampled))),
  rep("dementia", length(Cells(seurat_56_downsampled))),
  rep("no dementia", length(Cells(seurat_57_downsampled))),
  rep("dementia", length(Cells(seurat_58_downsampled))),
  rep("no dementia", length(Cells(seurat_59_downsampled))),
  rep("dementia", length(Cells(seurat_60_downsampled))),
  rep("no dementia", length(Cells(seurat_61_downsampled))),
  rep("no dementia", length(Cells(seurat_62_downsampled))),
  rep("no dementia", length(Cells(seurat_63_downsampled))),
  rep("dementia", length(Cells(seurat_64_downsampled))),
  rep("dementia", length(Cells(seurat_65_downsampled))),
  rep("no dementia", length(Cells(seurat_66_downsampled))),
  rep("dementia", length(Cells(seurat_67_downsampled))),
  rep("dementia", length(Cells(seurat_68_downsampled))),
  rep("no dementia", length(Cells(seurat_69_downsampled))),
  rep("dementia", length(Cells(seurat_70_downsampled))),
  rep("dementia", length(Cells(seurat_71_downsampled))),
  rep("dementia", length(Cells(seurat_72_downsampled))),
  rep("dementia", length(Cells(seurat_73_downsampled))),
  rep("no dementia", length(Cells(seurat_74_downsampled))),
  rep("no dementia", length(Cells(seurat_75_downsampled))),
  rep("dementia", length(Cells(seurat_76_downsampled))),
  rep("no dementia", length(Cells(seurat_77_downsampled))),
  rep("dementia", length(Cells(seurat_78_downsampled))),
  rep("dementia", length(Cells(seurat_79_downsampled))),
  rep("dementia", length(Cells(seurat_81_downsampled))),
  rep("dementia", length(Cells(seurat_82_downsampled))),
  rep("dementia", length(Cells(seurat_83_downsampled))),
  rep("no dementia", length(Cells(seurat_84_downsampled))),
  rep("no dementia", length(Cells(seurat_85_downsampled))),
  rep("dementia", length(Cells(seurat_86_downsampled))),
  rep("dementia", length(Cells(seurat_87_downsampled))),
  rep("dementia", length(Cells(seurat_88_downsampled))),
  rep("dementia", length(Cells(seurat_89_downsampled))),
  rep("dementia", length(Cells(seurat_90_downsampled))),
  rep("no dementia", length(Cells(seurat_91_downsampled)))
)

combined_data$condition <- conditions

head(combined_data$condition)
table(combined_data$condition)
Idents(combined_data) <- "condition" 
table(combined_data$condition)

Layers(combined_data)

combined_data <- JoinLayers(combined_data)
de_genes <- FindMarkers(combined_data, 
                        ident.1 = "dementia",  
                        ident.2 = "no dementia",  
                        test.use = "wilcox",  
                        logfc.threshold = 0.25, 
                        min.pct = 0.1)

significant_genes <- de_genes[de_genes$p_val_adj < 0.05 & abs(de_genes$avg_log2FC) > 0.1, ]
#0.1

print(significant_genes)
library(EnhancedVolcano)

EnhancedVolcano(de_genes, 
                lab = rownames(de_genes), 
                x = 'avg_log2FC', 
                y = 'p_val_adj')


top_genes <- rownames(significant_genes)[1:20]

print(top_genes)
DoHeatmap(combined_data, features = top_genes, group.by = "condition")
