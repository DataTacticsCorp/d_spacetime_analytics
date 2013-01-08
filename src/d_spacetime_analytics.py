#!/usr/bin/python

from utils import d_spacetime_analytics_utils
import numpy as np
import pysal
import os

class d_spacetime_analytics:
      
    def __init__(self, config):
        self.d_spacetime_analytics_main(config)
      
    def d_spacetime_analytics_main(self, config):
        #initialize utils class 
        self.utils = d_spacetime_analytics_utils()
    
        os.chdir(config['data_dir'])        #change the working directory to the configured data directory
        f = open(config['input_csv'],'r')   #open the input data csv
        lines = f.readlines()               #read all the lines into a list
        f.close()                           #close the file

        #Strip any quotes and split each line- each line is a list
        #Example data:
        #    ID,lat,lon,COUNTt1,COUNTt2
        #    11,3.58187154534,7.4475086296,4,0
        #    28,3.58187154534,7.44943566461,3,0
        #    32,3.58187154534,7.4498890846,12,0
        #    34,3.58187154534,7.4501157946,9,0
        lines = [line.strip().replace("\"", "").split(",") for line in lines]
        #create arrays for each field
        ids = [line[0] for line in lines[1:]]
        lat = [line[1] for line in lines[1:]]
        lon = [line[2] for line in lines[1:]]
        #data is a numpy array of the counts (last 2 columns)
        data = np.array([map(float,line[3:]) for line in lines[1:]])
    
        shp_data = {}
        shp_data['column_names'] = ['id', 'countt1', 'countt2']
        shp_data['column_types'] = ['N', 'N', 'N']
        for x in range(len(ids)):
            shp_data[ids[x]] = [lat[x], lon[x], int(data[x][0]), int(data[x][1])]

        self.utils.write_point_shapefile(config, shp_data)
        config['input_weights'] = self.utils.create_weights_from_shp(config)
        self.utils.fix_weights_file_header(config)
        
        #t1 and t2 are sorted lists for each count column, both are sent to the utils class to calculate their averages
        #The averages are then added to a numpy array and used to normalize the count columns resulting in Y
        t1 = sorted([float(line[3]) for line in lines[1:]])
        t2 = sorted([float(line[4]) for line in lines[1:]])
        adv = np.array([self.utils.get_mean_from_array(t1), self.utils.get_mean_from_array(t2)])
        Y = data/(adv*1.)
    
        #Open the weights file (gwt) using pysal- this requires the SHP/DBF files
        gwt = pysal.open(config['input_weights'])
        w = gwt.read()

        #pass the normalized (Y) numpy array, the weights file array (w) and set/pass the k value (number of circular sectors) to the rose function in utils
        #cuts is a list of values (in radians) defining circular sectors
        #theta, dx, dy  are arrays
        cuts, theta, dx, dy = self.utils.rose(Y, w, k=8)
        
        #z score calculation for both dx and dy- returns a numpy array
        zdx = self.utils.get_z_score_from_array(dx)
        zdy = self.utils.get_z_score_from_array(dy)
        
        #calculates the difference between zdx and zdy, returns a numpy array that is converted to a list
        #make a sorted copy in sorted_abs_dist
        abs_dist = self.utils.sqrt_list(np.ndarray.tolist(self.utils.get_abs_distance_of_zdx_zdy(zdx, zdy)))
        sorted_abs_dist = sorted(abs_dist)
        
        #calculate the following: min, max, median (50 percentile, Q2), 25 percentile (Q1), 75 percentile (Q3)
        percentile_0  = self.utils.get_percentile(sorted_abs_dist, 0.0, index=False)
        percentile_25 = self.utils.get_percentile(sorted_abs_dist, 0.25, index=False)
        percentile_50 = self.utils.get_percentile(sorted_abs_dist, 0.5, index=False)
        percentile_75 = self.utils.get_percentile(sorted_abs_dist, 0.75, index=False)
        percentile_100= self.utils.get_percentile(sorted_abs_dist, 1.0, index=False)
        #IQR = Q3 - Q1
        IQR = percentile_75 - percentile_25
        IQR_discontinuity = IQR * 1.5 + percentile_75
        IQR_EXdiscontinuity = IQR * 3.0 + percentile_75
        
        #create a dictionary for categories, define the categories with an ID and a range, append categories to itself
        categories = {}
        categories = self.utils.build_range_categories(categories, '0', [percentile_0, percentile_75])
        categories = self.utils.build_range_categories(categories, '1', [percentile_75, IQR_discontinuity])
        categories = self.utils.build_range_categories(categories, '2', [IQR_discontinuity, IQR_EXdiscontinuity])
        categories = self.utils.build_range_categories(categories, '3', [IQR_EXdiscontinuity, (percentile_100)])
        
        #sends variables to be printed to stdout
        self.utils.print_summary(IQR, IQR_discontinuity, IQR_EXdiscontinuity, categories, percentile_0, percentile_25, percentile_50, percentile_75, percentile_100)
    
        #open the file in the config dictionary to writing
        output = open(config['output_csv'], 'w')
        #write the header, then loop through to write the remainder of the data to CSV
        output.write("ids,lat,lon,t1_count,t2_count,theta,dx,dy,zdx,zdy,abs_dist,percentile_cat,lisa_cat\n")
        for x in range(len(theta)):
            #determine the category theta and abs_dist fall into
            lisa_cat = self.utils.get_category_from_range(self.utils.determin_lisa_category(cuts), theta[x])
            percentile_cat = self.utils.get_category_from_range(categories, abs_dist[x])
            output.write(str(ids[x])+","+str(lat[x])+","+str(lon[x])+","+str(t1[x])+","+str(t2[x])+","+str(theta[x])+","+str(dx[x])+","+str(dy[x])+","+str(zdx[x])+","+str(zdy[x])+","+str(abs_dist[x])+","+str(percentile_cat)+","+str(lisa_cat)+"\n")
        output.close()
        
if __name__ == "__main__": 
    #Configuration dictionary
    config = {}
    config['data_dir'] = "/home/ted/Desktop/Dev/workspace/d_spacetime_analytics/data"   #Directory where input data (csv, shp, dbf and gwt files) are located
    config['input_csv'] = "sample.csv"                                                #Input CSV
    config['output_csv'] = "sample_output.csv"                                        #Output file name to be written to
    config['output_shp'] = "sample_shp"
    config['shapefile_id_name'] = 'id'
    config['number_of_neighbors'] = 8
    config['input_weights'] = 'null'                                                    #Input GWT, set later

    d_spacetime_analytics(config)

