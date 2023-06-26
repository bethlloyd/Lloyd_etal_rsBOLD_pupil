
# Step 1) Make a template from the T1s. Mind you this is a fidily step. Tried multiple times without succes. Lastly, cd into T1 path. Did not make a template folder (ANTs created it). And moved the "bin" folder from the opt/AANTs location into the T1 path. Lastly, altered the "antsMultivariateTemplateConstruction2.sh" to make it work via qsub (in the end it took around ~24h) by removing "-q preemt" at 2 locations and increasing walltime to 72h and mem to 32gb.

inputPath=${PWD}/
outputPath=${PWD}/TemplateMultivariateBSplineSyN/

/home/memory/lindvoo/ANTs_NYU_LC/T1_im_orig/bin/antsMultivariateTemplateConstruction2.sh \
  -d 3 \
  -o ${outputPath}T_ \
  -i 4 \
  -g 0.2 \
  -j 2 \
  -c 4 \
  -k 1 \
  -w 1 \
  -f 8x4x2x1 \
  -s 3x2x1x0 \
  -q 100x70x50x10 \
  -n 1 \
  -r 1 \
  -l 1 \
  -m CC[2] \
  -t SyN \
  ${inputPath}/MRI*.nii

# Step 2) Delete everything ANTs created, EXCEPT the "T_template0.nii.gz" file which is the amazing template. Now we will register each T1 to this template (mind you within the first step this was also done but these are not optimale).

inputPath=${PWD}/
outputPath=${PWD}/TemplateMultivariateBSplineSyN/

/home/memory/lindvoo/ANTs_NYU_LC/T1_im_orig/bin/antsRegistrationSyN.sh \
  -d 3 \
  -f ${outputPath}T_template0.nii.gz \
  -m MRI_FCWML001_0001.nii \
  -o warped_



