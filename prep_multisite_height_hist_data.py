


bdir = '/Carnegie/DGE/Data/Shared/Labs/Asner/Private/Research/Researcher/Davies/4.Kruger/2.Termites_and_landuse/'

mound_center_csvs = [\
os.path.join(bdir,'shapefiles','andover_welv.csv'),\
os.path.join(bdir,'shapefiles','erosionOliver.csv'),\
os.path.join(bdir,'shapefiles','L1.csv'),\
os.path.join(bdir,'shapefiles','L23.csv'),\
os.path.join(bdir,'shapefiles','L45.csv'),\
os.path.join(bdir,'shapefiles','L7.csv'),\
os.path.join(bdir,'shapefiles','L8.csv'),\
os.path.join(bdir,'shapefiles','nwaswitshaka.csv')]

polygon_files = [\
os.path.join(bdir,'polygons','andover_raster.tif'),\
os.path.join(bdir,'polygons','erosion_raster.tif'),\
os.path.join(bdir,'polygons','landuse_1_raster.tif'),\
os.path.join(bdir,'polygons','landuse_2_raster.tif'),\
os.path.join(bdir,'polygons','landuse_456_raster.tif'),\
os.path.join(bdir,'polygons','landuse_7_raster.tif'),\
os.path.join(bdir,'polygons','landuse_8_raster.tif'),\
os.path.join(bdir,'polygons','nwas_raster.tif')]



key_df = pd.read_csv(bdir + '/polygons/polygon_key.csv',sep=',')
poly_size_numbers = []
poly_size_area = []

un_treat = np.unique(np.array(key_df['treatment']))
un_treat_label = ['Subsistance Agriculture','Communal Grazing', 'Kruger NP','Private Reserve']
height_list = [np.zeros(0) for x in range(len(un_treat))]
for _f in range(0,len(ensemble_files)):
 print(('filepair',mound_center_csvs[_f],polygon_files[_f]))

 df = pd.read_csv(mound_center_csvs[_f],sep=',')

 poly_set = gdal.Open(polygon_files[_f],gdal.GA_ReadOnly)
 poly = poly_set.ReadAsArray()

 un_poly = np.unique(poly)
 un_poly = un_poly[un_poly != 0]
 un_poly = un_poly[un_poly != 31]
 for m in range(0,len(un_poly)):
   poly_size_numbers.append(un_poly[m])
   poly_size_area.append(np.sum(poly == un_poly[m]))


 trans = poly_set.GetGeoTransform()
 
 per_ha_conv = 1.0e4/float(trans[1]*trans[1])
 
 px_y = ((np.array(df['y'])-trans[3])/float(trans[5])).astype(int)
 px_x = ((np.array(df['x'])-trans[0])/float(trans[1])).astype(int)

 #un_poly = np.unique(poly)
 #un_poly = un_poly[un_poly != 0]
 #for m in range(0,len(un_poly)):
 #  current_poly = un_poly[m]
 #  ret_list = get_mound_heights(poly == un_poly[m],df,px_x,px_y)
 #  height_list[un_treat.tolist().index(np.array(key_df['treatment'])[np.where(key_df['polygon_number'] == un_poly[m])][0])] = np.append(height_list[un_treat.tolist().index(np.array(key_df['treatment'])[np.where(key_df['polygon_number'] == un_poly[m])][0])],ret_list)

 #  np.savez('munged_dat/height_poly_' + str(un_poly[m]) + '.npz',histogram=ret_list)

np.savez('munged_dat/poly_areas.npz',poly_size_area = poly_size_area,poly_size_numbers=poly_size_numbers)
np.save('munged_dat/height_values.npy',height_list)
