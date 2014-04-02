cd IASP_50s_fwd/
cd MZZ
    echo 'in MZZ'
    ../../UTILS/xfield_transform > OUTPUT_field_transform 
cd ..

cd MXX_P_MYY/
    echo 'in MXX_P_MYY'
    ../../UTILS/xfield_transform > OUTPUT_field_transform  
cd ..

cd MXY_MXX_M_MYY/
    echo 'in MXY_MXX_M_MYY'
    ../../UTILS/xfield_transform > OUTPUT_field_transform  
cd ..

cd MXZ_MYZ/
    echo 'in MXZ_MYZ'
    ../../UTILS/xfield_transform > OUTPUT_field_transform  
cd ..

cd ../IASP_50s_bwd/
    echo 'in backward'
    ../UTILS/xfield_transform > OUTPUT_field_transform & 
cd ..

