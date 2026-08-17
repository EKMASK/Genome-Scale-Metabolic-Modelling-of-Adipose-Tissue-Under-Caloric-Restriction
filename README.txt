% 8CM00 Systems Medicine — Module 1, Group Assignment 3
% Project Structure
%
% 8CM00_Assignment3/
% ├── main_report.tex          ← Main LaTeX document (compile this)
% ├── README.txt               ← This file
% │
% ├── figures/                  ← MATLAB-generated figures go here
% │   └── (empty — populated by running MATLAB scripts)
% │       Q4_model_sizes_recon3d.png
% │       Q4_jaccard_heatmap_recon3d.png
% │       Q4_clustered_jaccard_recon3d.png
% │       Q4_dendrogram_recon3d.png
% │       Q4_MDS_jaccard_recon3d.png
% │       Q5_cosine_heatmap_recon3d.png
% │       Q5_clustered_cosine_recon3d.png
% │       Q5_cosine_dendrogram_recon3d.png
% │       Q5_MDS_cosine_recon3d.png
% │       Q5_objective_values_recon3d.png
% │
% ├── code/                     ← All MATLAB scripts
% │   ├── setup_models.m        ← Shared config: loads models + metadata
% │   ├── Q1.m                          ← Q1: PCA analysis
% │   ├── Q2.m                          ← Q2: context-specific models from Recon3D      
% │   │                                       (baseline)
% │   ├── Q3.m                          ← Q3: context-specific models from Foguet       
% │   │                                       Adipocyte Reconstruction (baseline)
% │   ├── Q4_jaccard_comparison.m       ← Q4: Jaccard similarity
% │   ├── Q5_FBA_cosine_similarity.m    ← Q5: FBA + cosine similarity
% │   ├── Q6.m                          ← Q6: context-specific models from Foguet       
% │   │                                       Adipocyte Reconstruction (post-           
% │   │                                       intervention)
% │   └── generate_iMAT_models.m        ← Backup: generate iMAT models
% │
% └── data/                     ← Expression data, models, measurements
%     ├── Adipose_Expression_WholeGenome.txt
%     ├── Adipose_Expression_Recon3.txt
%     ├── Adipose_Expression_Fouget.txt
%     ├── Yoyo_additional_measurements.txt
%     ├── Recon3DModel.mat
%     ├── Foguet_adipocyte.mat
%     ├── models_recon3d/       ← iMAT models from recon3d, unzip the file inside and   
%     │   │                       place them as follow
%     │   └── model_P01_baseline.mat
%     │   └── model_P02_baseline.mat, ...
%     └── models_fouget/        ← iMAT models from Fouget, unzip the files inside and   
%         │                       place them as follow (note that there are two zip     
%         │                       files inside)
%         └── model_P01_baseline.mat, 
%         └── model_P02_baseline.mat, ...
%
%
% HOW TO RUN
% ----------
% 1. Unzip the zip files inside models_fouget and recon3d files and place them as       
%    described above. You will put in total of 20 baseline models for recon3d in        
%    models_recon3d folder and 20 baseline plus 20 post-intervention models in          
%    models_fouget folder.
%    Place them in data/models_recon3d/ and/or data/models_fouget/.
%    Expected naming: model_P01_baseline.mat, model_P02_baseline.mat, etc.
%
% 2. Open MATLAB, cd to the code/ directory:
%       >> cd('/path/to/8CM00_Assignment3/code')
%
% 3. In code/setup_models.m, set reconstructionType to 'recon3d' or 'fouget'.
% 
% 4. Run Q1:
%       >> Q1
%
% 5. Run Q2:
%       >> Q2
%    Context-specific models using Recon3D for baseline samples will be saved in        
%    /data/iMAT_recon3_baseline_models/. The folder will be automatically created.
%
% 6. Run Q3:
%       >> Q3
%    Context-specific models using Foguet for baseline samples will be saved in         
%    /data/iMAT_foguet_baseline_models/. The folder will be automatically created.
% 
% 7. Run Q4:
%       >> Q4_jaccard_comparison
%    Figures are saved automatically to figures/.
%
% 8. Run Q5:
%       >> Q5_FBA_cosine_similarity
%    Figures are saved automatically to figures/.
%
% 9. Run Q6:
%       >> Q6
%    Context-specific models using Foguet for post-intervention samples will be saved in
%    /data/iMAT_foguet_post_models/. The folder will be automatically created.
%
% 10. Repeat steps 3-5 with the other reconstructionType to compare both.
%
%