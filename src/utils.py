#!/usr/bin/python

import numpy as np
import shapefile
import pysal as ps
import math
import os

class d_spacetime_analytics_utils:
        
    def debug_to_file(self, filename, data_list):
        '''
        Used for debugging any type of list
        @param filename: name of the file to write to
        @param data_list: list to be dumped to the file  
        '''
        data = ""
        for x in data_list:
            data += str(x)+"\n"
        out = open(filename, 'w')
        out.write(data)
        out.close()
        
    def get_mean_from_array(self, a):
        '''
        Simple function to calculate the mean from a list
        @parameter a: 1D array/list of values
        @return float of the mean
        '''
        size = len(a)
        total = 0.0
        for x in a:
            total += float(x)
        return (total/size)

    def sqrt_list(self, a):
        '''
        Square root all the values in a list
        @param a: any list/array
        @return: index preserved list of the sqrt of each entry in a
        '''
        sqrt_list = []
        for x in range(len(a)):
            sqrt_list.append(math.sqrt(a[x]))
        return sqrt_list

    def build_range_categories(self, categories, category_name, range_list):
        '''
        Build a dictionary of categories where the Key is the category name and it's value is the range list
        @param categories: a dictionary object (could be empty) to append the new category and range too
        @param category_name: name of the new category (Key)
        @param range_list: a list of ranges to set as the Value
        @return: categories dictionary
        '''
        categories[category_name] = range_list
        return categories

    def get_percentile(self, a, percent, index):
        '''
        if index = True the list index is returned that contains the percentile value
        if index = False the percentile value is calculated from the list of values
        '''
        if index:
            return self.get_percentile_index(a, percent)
        else:
            return self.get_percentile_value(sorted(a), percent)
        
    def get_percentile_index(self, a, percent):
        lenght = len(a)
        return lenght*percent

    def get_percentile_value(self, a, percent, key=lambda x:x):
        """
        Find the percentile of a list of values.
        @parameter a - is a list of values. Note N MUST BE already sorted.
        @parameter percent - a float value from 0.0 to 1.0.
        @parameter key - optional key function to compute value from each element of N.
        @return - the percentile of the values
        """
        if not a:
            return None
        k = (len(a)-1) * percent
        f = math.floor(k)
        c = math.ceil(k)
        if f == c:
            return key(a[int(k)])
        d0 = key(a[int(f)]) * (c-k)
        d1 = key(a[int(c)]) * (k-f)
        return d0+d1

    
    def determin_lisa_category(self, in_cuts):
        """
        Setup the lisa categories based on the number of cuts. Cuts should only be 4 or 8
        @parameter in_cuts - cut ranges, only used for it's length
        @return - a dictionary of categories
        """
        number_of_cuts = len(in_cuts)-1
        if number_of_cuts == 4:
            categories = {'HH': [0.0, 1.57079633], 'LH': [1.57079633, 3.14159265],\
                          'LL': [3.14159265, 4.71238898], 'HL': [4.71238898, 6.28318531]}
            return categories
            
        elif number_of_cuts == 8:
            categories = {'Hh': [0.0, 0.78539816], 'hH': [0.78539816, 1.57079633],\
                          'lH': [1.57079633, 2.35619449], 'Lh': [2.35619449, 3.14159265],\
                          'Ll': [3.14159265, 3.92699082], 'lL': [3.92699082, 4.71238898],\
                          'hL': [4.71238898, 5.49778714], 'Hl': [5.49778714, 6.28318531]}
            return categories
        else:
            return "EQUAL_TO_CUT"
    
    def get_category_from_range(self, categories, value):
        '''
        returns the category the value falls within given a dictionary. ex: category 'X': [0.0, 10.0])
        if value is 5 'X' is returned
        '''
        for cat in categories:
            if value >= categories[cat][0] and categories[cat][1] >= value:
                return cat

    def write_point_shapefile(self, config, shp_data):
        os.chdir(config['data_dir'])
        
        w = shapefile.Writer(shapefile.POINT)
        
        column_names = shp_data['column_names']
        column_types = shp_data['column_types']
        
        shp_data.pop('column_names')
        shp_data.pop('column_types')

        for name in column_names:
            for type in column_types:
                w.field(name, type)
                break;

        for point in shp_data:
            w.point(float(shp_data[point][1]), float(shp_data[point][0]))
            w.record(point, shp_data[point][2], shp_data[point][3])
            
        w.save(config['output_shp'])
        
    def create_weights_from_shp(self, config):
        os.chdir(config['data_dir'])
        shapefile_name = config['output_shp']+".shp"
        w = ps.knnW_from_shapefile(shapefile_name, k=config['number_of_neighbors'], idVariable=config['shapefile_id_name'])
        w.transform = 'r'
        fo = ps.open(config['output_shp']+".gwt", "w")
        fo.write(w)
        fo.close()
        return config['output_shp']+".gwt"

    def fix_weights_file_header(self, config):
        id = config['shapefile_id_name']
        file_name = config['output_shp']+'.gwt'
        os.chdir(config['data_dir'])
        f = open(file_name, 'r')
        lines = f.readlines()
        f.close()
        line = lines[0].strip().split()
        if len(line) == 4:
            line[2] = config['output_shp']
            line[3] = id
        else:
            line[2] = config['output_shp']
            line.append(id)
        lines[0] = str(line[0])+"\t"+str(line[1])+"\t"+line[2]+"\t"+line[3]+"\n"
        f = open(file_name, 'w')
        for line in lines:
            f.write(line)
        f.close()

    #simple z score calculation- input a list (a) and returns a z score for each item in a numpy array
    def get_z_score_from_array(self, a):
        mu = np.mean(a,None)
        sigma = self.samplestd(a)
        return (np.array(a)-mu)/sigma
    
    def get_abs_distance_of_zdx_zdy(self, zdx, zdy):
        dif = (zdx-zdy)
        return (dif*dif)

    #prints a summary of the results
    def print_summary(self, IQR, IQR_discontinuity, IQR_EXdiscontinuity, categories, percentile_0, percentile_25, percentile_50, percentile_75, percentile_100):
        print "IQR: "+str(IQR)
        print "UpperHinge(1.5):  "+str(IQR_discontinuity)
        print "EUpperHinge(3.0): "+str(IQR_EXdiscontinuity)
        print "Ranges:"
        print "   0: "+str(categories['0'][0])+" -> "+str(categories['0'][1])
        print "   1: "+str(categories['1'][0])+" -> "+str(categories['1'][1])
        print "   2: "+str(categories['2'][0])+" -> "+str(categories['2'][1])
        print "   3: "+str(categories['3'][0])+" -> "+str(categories['3'][1])+"+"
        print "Percentiles:"
        print "   MIN: "+str(percentile_0)
        print "   Q1: "+str(percentile_25)
        print "   Q2: "+str(percentile_50)
        print "   Q3: "+str(percentile_75)
        print "   MAX: "+str(percentile_100)

    #modified (simplified) rose function from pysal.spatial_dynamics.directional.rose
    def rose(self, Y, w, k=8):
        sw = 2*np.pi/k 
        cuts = np.arange(0.0,2*np.pi+sw,sw)
        wY = ps.lag_spatial(w,Y)
        dx = Y[:,-1]-Y[:,0]
        dy = wY[:,-1]-wY[:,0]
        theta = np.arctan2(dy,dx)
        neg = theta < 0.0
        utheta = theta*(1-neg) + neg * (2*np.pi+theta)         
        return cuts, utheta, dx, dy
    
    def samplestd(self, a, axis=0):
        return np.sqrt(self.samplevar(a,axis))
    
    def samplevar(self, a, axis=0):
        a, axis = self._chk_asarray(a, axis)
        mn = np.expand_dims(np.mean(a, axis), axis)
        deviations = a - mn
        n = a.shape[axis]
        svar = self.ss(deviations,axis) / float(n)
        return svar
    
    def ss(self, a, axis=0):
        a, axis = self._chk_asarray(a, axis)
        return np.sum(a*a, axis)
    
    def _chk_asarray(self, a, axis):
        if axis is None:
            a = np.ravel(a)
            outaxis = 0
        else:
            a = np.asarray(a)
            outaxis = axis
        return a, outaxis
            
    def get_min_max_from_list(self, in_list):
        return min(in_list), max(in_list)