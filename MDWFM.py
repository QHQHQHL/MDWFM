def TifftoArray(tifFile):
    dataset = gdal.Open(tifFile) 
    #des = dataset.GetDescription() 
    #x_size = dataset.RasterXSize
    #y_size = dataset.RasterYSize
    #geotransform = dataset.GetGeoTransform()
    #projection = dataset.GetProjection()
    #driver = dataset.GetDriver()
    tif_array = np.array(dataset.ReadAsArray())
    #row, col = tif_array.shape
    arraytif = np.nan_to_num(tif_array)
    r1 = arraytif.astype("float")
    return r1
def WriteToRasterFile(out_array, out_name, driver, Col, Row, Bandcount, DataType, NoDataValue, Geotransform, Geoprojection):
    output_raster = driver.Create(out_name, Col, Row, DataType, gdal.GDT_Float64) 
        
    outBand = output_raster.GetRasterBand(1)
    # write the data
    outBand.WriteArray(out_array, 0, 0)
    
    # flush data to disk, set the NoData value and calculate stats
    outBand.FlushCache()
    outBand.SetNoDataValue(NoDataValue)
    output_raster.SetGeoTransform(Geotransform)
    output_raster.SetProjection(Geoprojection)
    return output_raster
def Correlation_coefficient(a,b):
    if len(a) != 0:
        a_avg = sum(a)/len(a)
        b_avg = sum(b)/len(b)
        #print(sum(a),len(a))
        cov_ab = sum([(x - a_avg)*(y - b_avg) for x,y in zip(a, b)])
        sq = math.sqrt(sum([(x - a_avg)**2 for x in a])*sum([(x - b_avg)**2 for x in b]))
        
        if sq != 0:
            corr_factor = ((cov_ab/sq)+1.0)/2.0
        else:
            corr_factor = 1.0
        return corr_factor
    else:
        corr_factor = 1.0
        return corr_factor

def PCR6(S_mass,T_mass,P_mass,Q_mass,f):
    water = np.array([S_mass,T_mass,P_mass,Q_mass])
    land = np.array([(1.0-S_mass),(1.0-T_mass),(1.0-P_mass),(1.0-Q_mass)])
    PCR6_f1 = water[0]*water[1]*water[2]*water[3]
    #PCR6_m1
    PCR6_m1_lll = (land[1]*land[2]*land[3])/((water[0]**f)+(land[1]**f)+(land[2]**f)+(land[3]**f))
    PCR6_m1_llw = (land[1]*land[2]*water[3])/((water[0]**f)+(land[1]**f)+(land[2]**f)+(water[3]**f))
    PCR6_m1_lwl = (land[1]*water[2]*land[3])/((water[0]**f)+(land[1]**f)+(water[2]**f)+(land[3]**f))
    PCR6_m1_wll = (water[1]*land[2]*land[3])/((water[0]**f)+(water[1]**f)+(land[2]**f)+(land[3]**f))
    PCR6_m1_wwl = (water[1]*water[2]*land[3])/((water[0]**f)+(water[1]**f)+(water[2]**f)+(land[3]**f))
    PCR6_m1_wlw = (water[1]*land[2]*water[3])/((water[0]**f)+(water[1]**f)+(land[2]**f)+(water[3]**f))
    PCR6_m1_lww = (land[1]*water[2]*water[3])/((water[0]**f)+(land[1]**f)+(water[2]**f)+(water[3]**f))
    PCR_m1_sum = PCR6_m1_lll + PCR6_m1_llw + PCR6_m1_lwl + PCR6_m1_wll + PCR6_m1_wwl + PCR6_m1_wlw + PCR6_m1_lww
    PCR_m1 = (water[0]**f)*water[0]*PCR_m1_sum
    #PCR6_m2
    PCR6_m2_lll = (land[0]*land[2]*land[3])/((water[1]**f)+(land[0]**f)+(land[2]**f)+(land[3]**f))
    PCR6_m2_llw = (land[0]*land[2]*water[3])/((water[1]**f)+(land[0]**f)+(land[2]**f)+(water[3]**f))
    PCR6_m2_lwl = (land[0]*water[2]*land[3])/((water[1]**f)+(land[0]**f)+(water[2]**f)+(land[3]**f))
    PCR6_m2_wll = (water[0]*land[2]*land[3])/((water[1]**f)+(water[0]**f)+(land[2]**f)+(land[3]**f))
    PCR6_m2_wwl = (water[0]*water[2]*land[3])/((water[1]**f)+(water[0]**f)+(water[2]**f)+(land[3]**f))
    PCR6_m2_wlw = (water[0]*land[2]*water[3])/((water[1]**f)+(water[0]**f)+(land[2]**f)+(water[3]**f))
    PCR6_m2_lww = (land[0]*water[2]*water[3])/((water[1]**f)+(land[0]**f)+(water[2]**f)+(water[3]**f))
    PCR_m2_sum = PCR6_m2_lll + PCR6_m2_llw + PCR6_m2_lwl + PCR6_m2_wll + PCR6_m2_wwl + PCR6_m2_wlw + PCR6_m2_lww
    PCR_m2 = (water[1]**f)*water[1]*PCR_m2_sum
    #PCR6_m3
    PCR6_m3_lll = (land[0]*land[1]*land[3])/((water[2]**f)+(land[0]**f)+(land[1]**f)+(land[3]**f))
    PCR6_m3_llw = (land[0]*land[1]*water[3])/((water[2]**f)+(land[0]**f)+(land[1]**f)+(water[3]**f))
    PCR6_m3_lwl = (land[0]*water[1]*land[3])/((water[2]**f)+(land[0]**f)+(water[1]**f)+(land[3]**f))
    PCR6_m3_wll = (water[0]*land[1]*land[3])/((water[2]**f)+(water[0]**f)+(land[1]**f)+(land[3]**f))
    PCR6_m3_wwl = (water[0]*water[1]*land[3])/((water[2]**f)+(water[0]**f)+(water[1]**f)+(land[3]**f))
    PCR6_m3_wlw = (water[0]*land[1]*water[3])/((water[2]**f)+(water[0]**f)+(land[1]**f)+(water[3]**f))
    PCR6_m3_lww = (land[0]*water[1]*water[3])/((water[2]**f)+(land[0]**f)+(water[1]**f)+(water[3]**f))
    PCR_m3_sum = PCR6_m3_lll + PCR6_m3_llw + PCR6_m3_lwl + PCR6_m3_wll + PCR6_m3_wwl + PCR6_m3_wlw + PCR6_m3_lww
    PCR_m3 = (water[2]**f)*water[2]*PCR_m3_sum
    #PCR6_m4
    PCR6_m4_lll = (land[0]*land[1]*land[2])/((water[3]**f)+(land[0]**f)+(land[1]**f)+(land[2]**f))
    PCR6_m4_llw = (land[0]*land[1]*water[2])/((water[3]**f)+(land[0]**f)+(land[1]**f)+(water[2]**f))
    PCR6_m4_lwl = (land[0]*water[1]*land[2])/((water[3]**f)+(land[0]**f)+(water[1]**f)+(land[2]**f))
    PCR6_m4_wll = (water[0]*land[1]*land[2])/((water[3]**f)+(water[0]**f)+(land[1]**f)+(land[2]**f))
    PCR6_m4_wwl = (water[0]*water[1]*land[2])/((water[3]**f)+(water[0]**f)+(water[1]**f)+(land[2]**f))
    PCR6_m4_wlw = (water[0]*land[1]*water[2])/((water[3]**f)+(water[0]**f)+(land[1]**f)+(water[2]**f))
    PCR6_m4_lww = (land[0]*water[1]*water[2])/((water[3]**f)+(land[0]**f)+(water[1]**f)+(water[2]**f))
    PCR_m4_sum = PCR6_m4_lll + PCR6_m4_llw + PCR6_m4_lwl + PCR6_m4_wll + PCR6_m4_wwl + PCR6_m4_wlw + PCR6_m4_lww
    PCR_m4 = (water[3]**f)*water[3]*PCR_m4_sum
    PCR6 = PCR6_f1 + PCR_m1 + PCR_m2 + PCR_m3 + PCR_m4
    return PCR6
    def ImageFusion(ls_path,s1_path,s2_path,S_LS,S_S1,S_S2,J,f,th):
    #Remote-sensing images for fusion
    ls = TifftoArray(ls_path)
    s1 = TifftoArray(s1_path)
    s2 = TifftoArray(s2_path)
    
    #construction of the fusion result
    result = TifftoArray(s2_path)
    row, col = s2.shape
    print(row,col)
    for i in range(row):
        for j in range(col):
            result[i][j] = 0.0
    resultPCR6 = TifftoArray(s2_path)
    row, col = s2.shape
    print(row,col)
    for i in range(row):
        for j in range(col):
            resultPCR6[i][j] = 0.0
    #construction of the multiple dimension
    resultS = TifftoArray(s2_path)
    row, col = s2.shape
    for i in range(row):
        for j in range(col):
            resultS[i][j] = 0.0
    resultT = TifftoArray(s2_path)
    row, col = s2.shape
    for i in range(row):
        for j in range(col):
            resultT[i][j] = 0.0
    resultP = TifftoArray(s2_path)
    row, col = s2.shape
    for i in range(row):
        for j in range(col):
            resultP[i][j] = 0.0
    resultQ = TifftoArray(s2_path)
    row, col = s2.shape
    for i in range(row):
        for j in range(col):
            resultQ[i][j] = 0.0
    #Property value for each image (spatial resolution)
    pvs1 = float(1.0/10.0)**2.0
    pvs2 = float(1.0/10.0)**2.0
    pvls = float(1.0/30.0)**2.0
    #Quality value for each image (training sample separability)
    qvs1 = S_S1**2.0
    qvs2 = S_S2**2.0
    qvls = S_LS**2.0
    #Spatial value and temporal value for each pixel
    for i in range(J,row-J):
        for j in range(J,col-J):
            if s2[i][j] != 2.0 and ls[i][j] != 2.0:
                S_mass = Spatialvalue(s2, i, j, 1, 3.0) #water judgment of spatial dimension
                # obtain 3D neighboring pixels according to J
                ns2_ls = [s2[i][j]]
                ns2_s1 = [s2[i][j]]
                ns1 = [s1[i][j]]
                nls = [ls[i][j]]
                for nn in range(J):
                    ns2_ls = ns2_ls + Time_Neighbor_OPOP(s2, nn+1, i, j, s2, ls)
                    ns2_s1 = ns2_s1 + Time_Neighbor_SAROP(s2, nn+1, i, j, s2)
                    ns1 = ns1 + Time_Neighbor_SAROP(s1, nn+1, i, j, s2)
                    nls = nls + Time_Neighbor_OPOP(ls, nn+1, i, j, s2, ls)
                #Temporal value for the pixels with same position in different images
                tvs2 = 1.0**1.0
                tvs1 = Correlation_coefficient(ns2_s1, ns1)**1.0
                tvls = Correlation_coefficient(ns2_ls, nls)**1.0
                #Temporal weight
                T_weighted_s1 = float(tvs1/(tvs1+tvs2+tvls))
                T_weighted_s2 = float(tvs2/(tvs1+tvs2+tvls))
                T_weighted_ls = float(tvls/(tvs1+tvs2+tvls))
                T_mass = T_weighted_ls*ls[i][j]+(T_weighted_s1*s1[i][j])+(T_weighted_s2*s2[i][j]) #water judgment of temporal dimension
                #Property weight
                P_weighted_s1 = float(pvs1/(pvs1+pvs2+pvls))
                P_weighted_s2 = float(pvs2/(pvs1+pvs2+pvls))
                P_weighted_ls = float(pvls/(pvs1+pvs2+pvls))
                P_mass = P_weighted_ls*ls[i][j]+(P_weighted_s1*s1[i][j])+(P_weighted_s2*s2[i][j]) #water judgment of Property dimension
                #Quality weight
                Q_weighted_s1 = float(qvs1/(qvs1+qvs2+qvls))
                Q_weighted_s2 = float(qvs2/(qvs1+qvs2+qvls))
                Q_weighted_ls = float(qvls/(qvs1+qvs2+qvls))
                Q_mass = Q_weighted_ls*ls[i][j]+(Q_weighted_s1*s1[i][j])+(Q_weighted_s2*s2[i][j]) #water judgment of Quality dimension
                resultS[i][j] = S_mass#PCR6(S_mass,T_mass,P_mass,Q_mass,f)# image fusion
                resultT[i][j] = T_mass
                resultP[i][j] = P_mass
                resultQ[i][j] = Q_mass
                resultPCR6[i][j] = PCR6(S_mass,T_mass,P_mass,Q_mass,f)# image fusion
                result[i][j] = 1.0 - (S_mass*T_mass*P_mass*Q_mass) - ((1.0-S_mass)*(1.0-T_mass)*(1.0-P_mass)*(1.0-Q_mass)) # Conflict
                #print(S_mass,T_mass,P_mass,Q_mass,result[i][j]) # Conflict
            if s2[i][j] != 2.0 and ls[i][j] == 2.0:
                S_mass = Spatialvalue(s2, i, j, 1, 3.0) #water judgment of spatial dimension
                # obtain 3D neighboring pixels according to J
                ns2_s1 = [s2[i][j]]
                ns1 = [s1[i][j]]
                for nn in range(J):
                    ns2_s1 = ns2_s1 + Time_Neighbor_SAROP(s2, nn+1, i, j, s2)
                    ns1 = ns1 + Time_Neighbor_SAROP(s1, nn+1, i, j, s2)
                tvs2 = 1.0**1.0
                tvs1 = Correlation_coefficient(ns2_s1, ns1)**1.0
                #Temporal weight
                T_weighted_s1 = float(tvs1/(tvs1+tvs2))
                T_weighted_s2 = float(tvs2/(tvs1+tvs2))
                T_mass = (T_weighted_s1*s1[i][j])+(T_weighted_s2*s2[i][j]) #water judgment of temporal dimension
                #Property weight
                P_weighted_s1 = float(pvs1/(pvs1+pvs2))
                P_weighted_s2 = float(pvs2/(pvs1+pvs2))
                P_mass = (P_weighted_s1*s1[i][j])+(P_weighted_s2*s2[i][j]) #water judgment of Property dimension
                #Quality weight
                Q_weighted_s1 = float(qvs1/(qvs1+qvs2))
                Q_weighted_s2 = float(qvs2/(qvs1+qvs2))
                Q_mass = (Q_weighted_s1*s1[i][j])+(Q_weighted_s2*s2[i][j]) #water judgment of Quality dimension
                resultS[i][j] = S_mass#PCR6(S_mass,T_mass,P_mass,Q_mass,f)
                resultT[i][j] = T_mass
                resultP[i][j] = P_mass
                resultQ[i][j] = Q_mass
                resultPCR6[i][j] = PCR6(S_mass,T_mass,P_mass,Q_mass,f)# image fusion
                result[i][j] = 1.0 - (S_mass*T_mass*P_mass*Q_mass) - ((1.0-S_mass)*(1.0-T_mass)*(1.0-P_mass)*(1.0-Q_mass))
                #print(S_mass,T_mass,P_mass,Q_mass,result[i][j])
            if s2[i][j] == 2.0 and ls[i][j] != 2.0:
                S_mass = Spatialvalue(s2, i, j, 1, 3.0)
                tvs1 = 0.5
                tvls = 0.5
                #Temporal weight
                T_weighted_s1 = float(tvs1/(tvs1+tvls))
                T_weighted_ls = float(tvls/(tvs1+tvls))
                T_mass = T_weighted_ls*ls[i][j]+(T_weighted_s1*s1[i][j]) #water judgment of temporal dimension
                #Property weight
                P_weighted_s1 = float(pvs1/(pvs1+pvls))
                P_weighted_ls = float(pvls/(pvs1+pvls))
                P_mass = P_weighted_ls*ls[i][j]+(P_weighted_s1*s1[i][j]) #water judgment of Property dimension
                #Quality weight
                Q_weighted_s1 = float(qvs1/(qvs1+qvls))
                Q_weighted_ls = float(qvls/(qvs1+qvls))
                Q_mass = Q_weighted_ls*ls[i][j]+(Q_weighted_s1*s1[i][j]) #water judgment of Quality dimension
                if S_mass != 2.0:
                    resultS[i][j] = S_mass#PCR6(S_mass,T_mass,P_mass,Q_mass,f)
                    resultT[i][j] = T_mass
                    resultP[i][j] = P_mass
                    resultQ[i][j] = Q_mass
                    resultPCR6[i][j] = PCR6(S_mass,T_mass,P_mass,Q_mass,f)# image fusion
                    result[i][j] = 1.0 - (S_mass*T_mass*P_mass*Q_mass) - ((1.0-S_mass)*(1.0-T_mass)*(1.0-P_mass)*(1.0-Q_mass))
                    #print(S_mass,T_mass,P_mass,Q_mass,result[i][j])
                if S_mass == 2.0:
                    resultS[i][j] = 0.0#PCR6_3D(T_mass,P_mass,Q_mass,f)
                    resultT[i][j] = T_mass
                    resultP[i][j] = P_mass
                    resultQ[i][j] = Q_mass
                    resultPCR6[i][j] = PCR6_3D(T_mass,P_mass,Q_mass,f)# image fusion
                    result[i][j] = 1.0 - (T_mass*P_mass*Q_mass) - ((1.0-T_mass)*(1.0-P_mass)*(1.0-Q_mass))
                    #print(T_mass,P_mass,Q_mass,result[i][j])
            if ls[i][j] == 2.0 and s2[i][j] == 2.0:
                resultS[i][j] = 0.0#s1[i][j]
                resultT[i][j] = 0.0
                resultP[i][j] = 0.0
                resultQ[i][j] = 0.0
                result[i][j] = 0.0
                resultPCR6[i][j] = s1[i][j]# image fusion
    dataset = gdal.Open(s2_path)
    geotransform = dataset.GetGeoTransform()
    geoprojection = dataset.GetProjection()
    driver = dataset.GetDriver()
    array = np.array(dataset.ReadAsArray())
    row, col = array.shape
    dataType = 1
    noDataValue = 0
    bandcount = 1
    out_nameS = "E:/E/F1/suri/S_imagefusion"+"_f="+str(f)+"_"+str(th)+".tif"
    out_nameT = "E:/E/F1/suri/T_imagefusion"+"_f="+str(f)+"_"+str(th)+".tif"
    out_nameP = "E:/E/F1/suri/P_imagefusion"+"_f="+str(f)+"_"+str(th)+".tif"
    out_nameQ = "E:/E/F1/suri/Q_imagefusion"+"_f="+str(f)+"_"+str(th)+".tif"
    out_nameConf = "E:/E/F1/suri/Conf_imagefusion"+"_f="+str(f)+"_"+str(th)+".tif"
    out_namePCR6 = "E:/E/F1/suri/PCR6_imagefusion"+"_f="+str(f)+"_"+str(th)+".tif"
    outputrasterPCR6=WriteToRasterFile(resultPCR6, out_namePCR6, driver, col, row, bandcount, dataType, noDataValue, geotransform, geoprojection)
    outputraster=WriteToRasterFile(result, out_nameConf, driver, col, row, bandcount, dataType, noDataValue, geotransform, geoprojection)
    outputraster1=WriteToRasterFile(resultS, out_nameS, driver, col, row, bandcount, dataType, noDataValue, geotransform, geoprojection)
    outputraster2=WriteToRasterFile(resultT, out_nameT, driver, col, row, bandcount, dataType, noDataValue, geotransform, geoprojection)
    outputraster3=WriteToRasterFile(resultP, out_nameP, driver, col, row, bandcount, dataType, noDataValue, geotransform, geoprojection)
    outputraster4=WriteToRasterFile(resultQ, out_nameQ, driver, col, row, bandcount, dataType, noDataValue, geotransform, geoprojection)
    print(outputraster4)
    return(result)
