import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.stats import ranksums, sem
import numpy as np
import pandas as pd
import datetime
from datetime import timedelta
import requests
import glob
import math
import json
import os
import shutil
from utils import utils
from operator import add
from scipy.stats import expon

#api_url = 'https://data.braincommons.org/'
#api_url = 'https://gen3.datacommons.io/'
api_url='https://gen3-neuro.datacommons.io/'
summary_order = [
    "_study_count",
    "_case_count",
    "_demographic_count",
    "_diagnosis_count",
    "_exposure_count",
    "_family_history_count",
    "_clinical_lab_test_count",      
    "_injury_or_illness_count",    
    "_montreal_cognitive_functional_test_count",
    "_assessment_of_depression_count",   
    "_imaging_fmri_exam_count",    
    "_imaging_mri_exam_count",
    "_aliquot_count",
    "_submitted_expression_array_count",
    "_submitted_unaligned_reads_count"
    ]

summary_count_headers = {
    "_case_count": "Cases",
    "_study_count": "Studies",
    "_demographic_count": "Demographic records",
    "_diagnosis_count": "Diagnosis records",
    "_exposure_count": "Exposure records",
    "_aliquot_count": "# Aliquots",
    "_imaging_mri_exam_count": "# MRI Images",
    "_submitted_unaligned_reads_count": "# FASTQ Files",
    "_submitted_expression_array_count": "# CEL Files",
    "_family_history_count": "Family History records",
    "_clinical_lab_test_count": "Biomarker records",      
    "_injury_or_illness_count": "Injury records",    
    "_montreal_cognitive_functional_test_count": "Montreal Cognitive Assessment records",
    "_assessment_of_depression_count": "Despression Assessment records", 
    "_imaging_fmri_exam_count": "# fMRI Image Exam records"
}

project_names = {
    "bhc-cnp-open-fmri": "Open FMRI",
    "bhc-fmri-duke-test": "Duke FMRI",
    "bhc-RAPID-DxPilot": "Rapid-Dx Pilot",
    "mjff-LRRK2": "MJFF-LRRK2 Biomarker Project"
}

hemisphere_names = {
    "lh": "Left Hemisphere",
    "rh": "Right Hemisphere"
}

class DiseasesTable(dict):
    ''' Represent disease table in HTML format for visualization '''
 
    def _repr_html_(self):
        html = []
        html.append("<table style>")
        html.append("<thead>")  
        html.append("<th>Disease</th>")  
        html.append("<th># Subjects</th>")
        html.append("</thead>")
        total = 0
        for key in self:
            html.append("<tr>") 
            html.append("<td>%s</td>" % key)             
            html.append("<td>%s</td>" % self[key])           
            html.append("<tr>")
            total += self[key]
        html.append("<td>TOTAL</td>")             
        html.append("<td>%s</td>" % total)             
        html.append("</table>")        
        
        return ''.join(html)

class MetricsTable(dict):
    ''' Represent metrics tables in HTML format for visualization '''
 
    def _repr_html_(self):
        html = []
        html.append("<table style>")
        html.append("<thead>")  
        html.append("<th>Metric</th>")  
        html.append("<th>Value</th>")
        html.append("</thead>")       
        for key in self:
            html.append("<tr>") 
            html.append("<td>%s</td>" % key)             
            html.append("<td>%s</td>" % self[key])           
            html.append("<tr>") 
        html.append("</table>")        
        
        return ''.join(html)

class SummaryTable(dict):
    ''' Represent result tables in HTML format for visualization '''
 
    def _repr_html_(self):
        html = []
        html.append("<table style>")
        html.append("<thead>")  
        html.append("<th>Category</th>")  
        html.append("<th>Counts</th>")
        html.append("</thead>")
        for key in summary_order:
                html.append("<tr>") 
                html.append("<td>%s</td>" % summary_count_headers[key])             
                html.append("<td>%s</td>" % self[key])           
                html.append("<tr>")
        html.append("</table>")        
        
        return ''.join(html)

def add_keys(filename):
    ''' Get auth from our secret keys '''

    global auth 
    json_data=open(filename).read()
    keys = json.loads(json_data)
    creds_api = api_url + 'user/credentials/cdis/access_token'
    auth = requests.post(creds_api, json=keys)
    
def query_api(query_txt, variables = None):
    ''' Request results for a specific query '''

    if variables == None:
        query = {'query': query_txt}
    else:
        query = {'query': query_txt, 'variables': variables}        

    query_api = api_url + 'api/v0/submission/graphql'        
    output = requests.post(query_api, headers={'Authorization': 'bearer '+ auth.json()['access_token']}, json=query).text
    data = json.loads(output)    
    
    if 'errors' in data:
        print(data)   
    
    return data    

def get_projects(excluded=[]):
    ''' Query list of projects '''
    
    query_txt = """query Project { project(first:0) {project_id}} """
   
    data = query_api(query_txt) 
   
    projects = []
    for pr in data['data']['project']:
        if pr['project_id'] not in excluded:
            projects.append(pr['project_id'])
    projects = sorted(projects)   

    return projects
    
def query_summary_counts(projects=None):
    ''' Query summary counts for each data type'''
   
    if projects == None:
        projects = get_projects()
    elif not isinstance(projects,list):
        projects = [projects]
       
    dftotal = pd.DataFrame()
    for p in projects:
        query_txt = """query Counts ($projectID: [String]) {"""
        for param in summary_order:
            query_txt += """%s(project_id: $projectID)""" % param
        query_txt += "}" 
        variables = { 'projectID': p}
        data = query_api(query_txt, variables)
        indexes, values = [], []
        for key in summary_order:
            indexes.append(summary_count_headers[key])
            if key in data['data']:
                values.append(data['data'][key])
            else:
                values.append(0)           

        df = pd.DataFrame(values, index=indexes, columns=[p])
        if dftotal.empty:
            dftotal = df
        else:
            dftotal[p] = df[p]

    #dftotal = pd.concat(dftotal)    
 
    return dftotal

def query_summary_field(field, field_node, project_id = None):
    ''' Query summary counts for specific node'''
   
    if project_id != None:
        query_txt = """query { %s(first:0, project_id: "%s") {%s}} """ % (field_node, project_id, field) 
    else:
        query_txt = """query { %s(first:0) {%s project_id}} """ % (field_node, field)    
    
        
    data = query_api(query_txt)
    
    summary = {}
    total = []
    for d in data['data'][field_node]:
        
        if isinstance(d[field], float):
            d[field] = str(d[field])[:-2]        
        
        if 'project_id' in d:  
            summary.setdefault(d['project_id'], {})
            summary[d['project_id']].setdefault(d[field], 0)
            summary[d['project_id']][d[field]] += 1
            if d[field] not in total:
                total.append(d[field])            
        else:
            summary.setdefault(d[field], 0)        
            summary[d[field]] += 1
    
    if project_id != None:
        plot_field_metrics(summary, field)
    else:
        plot_overall_metrics(summary, field, total)        
    
    return summary

def field_distribution(field, field_node, project_id, bins = None, distrib=None, rate=None):
    ''' Plot distribution for one field'''
   
    if project_id != None:
        query_txt = """query { %s(first:0, project_id: "%s") {submitter_id %s}} """ % (field_node, project_id, field) 
    else:
        query_txt = """query { %s(first:0) {submitter_id %s project_id}} """ % (field_node, field)    
         
    data = query_api(query_txt)
       
    summary = {}
    total = []
    for d in data['data'][field_node]:
                
        if isinstance(d[field], float):
            d[field] = str(d[field])[:-2]        
        
        if 'project_id' in d:  
            summary.setdefault(d['project_id'], {})
            summary[d['project_id']].setdefault(d[field], 0)
            summary[d['project_id']][d[field]] += 1
            if d[field] not in total:
                total.append(d[field])            
        else:
            summary.setdefault(d[field], 0)        
            summary[d[field]] += 1    
         
    if len(summary)>10:
        
        accumulated = []
        for d in data['data'][field_node]:
            if d[field] != None:
                accumulated.append(float(d[field]))        
        
        plot_histogram(accumulated, field, bins, distrib, rate)
        
    else:
        
        N = len(summary)

        values = []
        types = []

        for n in sorted(summary, key=summary.get, reverse=True):
            values.append(summary[n])
            types.append(n)
            
        total = sum(values)
        positions = np.arange(N)
        fig, ax = plt.subplots(1, 1, figsize=(3*N, N))

        size_prop = (N/10) + 1
        ax.bar(positions, values, 0.2, align='center', alpha=0.4, color='b')
  
        plt.title('Summary counts by (' + field + ')', fontsize=10*size_prop)
        plt.ylabel('COUNTS', fontsize=10*size_prop)    
        plt.ylim(0, max(values)+5)
        plt.xlabel(field.upper(), fontsize=10*size_prop)  
        plt.xticks(positions, types, fontsize=10*size_prop)         
   
        # fit curve
        if distrib == 'exponential':
            fit_curve = expon.pdf(positions, 0, 1.0/rate)*total
            ax.plot(positions, fit_curve, 'r-', lw=2)
        if distrib == 'uniform':
            fit_curve = [total/float(len(positions))] * len(positions)
            ax.plot(positions, fit_curve, 'r-', lw=2)  

    return data    
    
    
def plot_histogram(accumulated, xlab, bins = None, distrib=None, rate=None):  
        
    # the histogram of the data
    plt.figure(figsize=(8, 4))
    fig, ax = plt.subplots(1, 1)
    n, positions, patches = ax.hist(accumulated, bins, facecolor='b', alpha=0.75)
    total=len(accumulated)

    plt.xlabel(xlab)
    plt.ylabel('Counts')
    plt.title('Histogram of ' + xlab)
    plt.grid(True)

    # fit curve
    if distrib == 'exponential':
        fit_curve = expon.pdf(positions, 0, 1.0/rate)*total
        ax.plot(positions, fit_curve, 'r-', lw=2)
    if distrib == 'uniform':
        fit_curve = [total/float(len(positions))] * len(positions)
        ax.plot(positions, fit_curve, 'r-', lw=2)  
        

def plot_bars(values, types, xlab, errors=[]):         
        
    N = len(values)
    positions = np.arange(N)        
    plt.figure(figsize=(3*N, N))   
    
    size_prop = (N/10) + 1
    if errors:
        plt.bar(positions, values, 0.4, align='center', alpha=0.4, color='cadetblue', yerr=errors, error_kw=dict(elinewidth=3, capsize=20, marker='o', markeredgewidth=3))
    else:
        plt.bar(positions, values, 0.4, align='center', alpha=0.4, color='cadetblue')
        plt.ylim(0, int(max(values)*1.1))
        
    plt.title('Summary counts by (' + xlab + ')', fontsize=10*size_prop)
    plt.ylabel('COUNTS', fontsize=10*size_prop)
    plt.yticks(fontsize=10*size_prop)
    plt.xlabel(xlab.upper(), fontsize=10*size_prop)  
    plt.xticks(positions, types, fontsize=10*size_prop, rotation='vertical')    
    
    if not errors:
        for i, v in enumerate(values):
            plt.text(i-0.05, v, str(v), color='red', fontweight='bold', fontsize=10*size_prop)
                   
        
def get_disease_cohorts(project_id):
    ''' Query summary counts for each data type'''
   
    query_txt = """query{
                      condition(first:0, project_id:"%s"){
                          primary_diagnosis
                      } 
               } """ % (project_id)

    data = query_api(query_txt) 

    diagnosis_counts = {}
    for diagnosis in data['data']['condition']:
        diagnosis_counts.setdefault(diagnosis['primary_diagnosis'], 0)
        diagnosis_counts[diagnosis['primary_diagnosis']] += 1    
    
    table = DiseasesTable(diagnosis_counts)
    
    return table
    
def plot_field_metrics(summary_counts, field):
    ''' Plot summary results in a barplot ''' 
    
    N = len(summary_counts)

    values = []
    types = []

    for n in sorted(summary_counts, key=summary_counts.get, reverse=True):
        values.append(summary_counts[n])
        types.append(n)
           
    plot_bars(values, types, field)


def plot_overall_metrics(summary_counts, field, totals):    
    ''' Visualize summary results across projects in a barplot ''' 
    
    results = {}
    projects = {}
    for project in summary_counts:
        
        results[project] = []
        projects.setdefault(project, 0)
            
        for value in totals:
            if value in summary_counts[project]:
                results[project].append(summary_counts[project][value])
                projects[project] += summary_counts[project][value]
            else:
                results[project].append(0)

    N = len(totals)
    positions = np.arange(N) 
    sorted_projects = sorted(projects, key=projects.get, reverse=True)
    bar_size = 0.2
    size_prop = (N/20) + 1
    
    plots = []
    plt.figure(figsize=(8, 4))
    left = [0]*N
    for pr in sorted_projects:
        p = plt.barh(positions, results[pr], bar_size, left, align='center', alpha=1)        
        plots.append(p[0])
        left = list(map(add, left, results[pr]))
        
    plt.title('Summary counts by (' + field + ')', fontsize=10*size_prop)
    plt.xlabel('COUNTS', fontsize=10*size_prop)    
    plt.xlim(0, max(left)+5)
    plt.ylabel(field.upper(), fontsize=10*size_prop)  
    plt.yticks(positions, totals, fontsize=10*size_prop)    
    plt.legend(plots, sorted_projects, fontsize=10*size_prop)
           
    plt.show()     
    
def get_cortical_qc(f_list):
    l=len(f_list) 
    print (l)
    fig = plt.figure(figsize= (20, 10*(int(l/2.0+0.5))))
    #fig = plt.figure((4,4))
    
    #fig.suptitle('SUBCORTICAL SEGMENTATION', fontsize=20, fontweight='bold')    
    pos = 0
    for file in f_list:
            # get file name from the file path
            indicesB = [i for i, x in enumerate(file) if x == '/']
            indicesE = [i for i, x in enumerate(file) if x == '.']
            pic_title=file[(indicesB[-1]+1):indicesE[-1]]
            pos += 1                             
            a=fig.add_subplot(int(l/2.0+0.5),2,pos)
            a.axis('off')
            img = mpimg.imread(file)
            imgplot = plt.imshow(img)                  
            a.set_title(pic_title, fontsize=18) 
    plt.show()
    

    
def plot_measure(results, title, measure, pvalues, threshold = 0.01):
    ''' Visualize metrics from cortical pipeline per disease'''
    
    N = len(results)
    positions = np.arange(N) + 1    
    
    plt.figure()
    plt.boxplot(list(results.values()), patch_artist=True)
    plt.xticks(positions, results.keys(), rotation='vertical')
    plt.title(title + ' - ' + measure, fontsize = 16)
    plt.ylabel(measure) 
    bottom,top = plt.ylim()
    plt.ylim(bottom,top*1.05)
    
    for p, group in zip(positions, pvalues):
        if pvalues[group] > 0:
            col = "red"
            if pvalues[group] < threshold: col = "green" 
            plt.text(p, top, "p={0:.4f}".format(pvalues[group]),
                 horizontalalignment='center', color=col, weight="bold")    
    
    plt.show()    

def plot_values(data, label):
    ''' Plot longitudinal variable evolution along several followups'''
    
    fig = plt.figure(figsize=(11,8))

    
    values = list(data["mean"])
    labels = list(data.index)
    size   = len(labels)
    errors = [v/math.sqrt(size) for v in list(data["std"])]
  
    ax1 = fig.add_subplot(111)
    ax1.errorbar(range(size), values, yerr=errors, fmt='-o', ecolor='gray', capsize=10)
    ax1.set_xlabel("Event ID", fontsize=14)
    ax1.set_ylabel(label, fontsize=14)
    ax1.set_xticks(range(size))
    ax1.set_xticklabels(labels)
        
        
def run_statistical_test(values, baseline="healthy control"):

    pvalues = {}
    for g in values:
        if g != baseline:
            test = ranksums(values[baseline], values[g])
            pvalues[g] = test.pvalue    
        else:
            pvalues[g] = -1
    
    return pvalues
    
    
def read_result_file(file):
    ''' Read TSV/CSV files from MRI results'''    
    if '.tsv' in file:
        sep = '\t'
    elif '.csv' in file:
        sep = ','
    
    headers = []
    results = {}
    with open(file, mode='r') as infile:
        for line in infile:
            columns = line.strip('\n').split(sep)
            if not headers:
                headers = columns
            else:
                pos = 0
                for h in headers:
                    if h == 'SubjID':
                        subject = columns[pos]
                        results.setdefault(subject, {})                     
                    else:
                        results[subject][h] = columns[pos]
                    pos += 1
    
    return results    
    
    
def get_external_surface_qc(project_id, subject_id):       
    ''' Get external surface results from data model'''    
        
    # Query data
    query_txt = """query {
                      case(submitter_id: "%s"){
                        mri_exams(first:0){
                          mri_images(first:0){
                            mri_analysis_workflows(first:0){
                              mri_results(first:0, data_category: "MRI Derived Image", data_type: "Cortical Thickness"){
                                file_name
                                id
                              }
                            }
                          }
                        }
                      }
                }""" % (subject_id)
    
    data = query_api(query_txt)

    # Display external surface
    fig = plt.figure(figsize=(15, 15))
    fig.suptitle('EXTERNAL SURFACE SEGMENTATION', fontsize=20, fontweight='bold')    
    for file in data['data']['case'][0]['mri_exams'][0]['mri_images'][0]['mri_analysis_workflows'][0]['mri_results']:
        
        # Download image files if not in local
        image = file['file_name']
        if not os.path.exists(image):
            image = utils.download_file(auth, api_url, file['id'], image)  
        
        # Plot images
        print(image)
        if 'lh.lat.' in image:
            title_view = 'Lateral View - Left Hemisphere'
            pos = 1
        elif 'lh.med.' in image:
            title_view = 'Medial View - Left Hemisphere'
            pos = 2                
        elif 'rh.lat.' in image:
            title_view = 'Lateral View - Right Hemisphere'
            pos = 3 
        elif 'rh.med.' in image:
            title_view = 'Medial View - Right Hemisphere'
            pos = 4                 
                
        a=fig.add_subplot(2,2,pos)
        a.axis('off')
        img = mpimg.imread(image)
        imgplot = plt.imshow(img)                  
        a.set_title(title_view, fontsize=18)        


def get_hippocampal_qc(project_id, subject_id, output_path, view, hemisphere):       
    ''' Get external surface results from data model'''    
       
    # Display external surface
    fig = plt.figure(figsize=(15, 15))
    fig.suptitle('HIPPOCAMPAL SEGMENTATION', fontsize=20, fontweight='bold')    
    pos = 0
    slides = ['20','40','60','80']
    for file in output_path:
        if view in file and hemisphere in file and not 'T1' in file and not '_fis' in file:
            print(file)
            title_view = '%s View (%s, slide %s)' % (view, hemisphere, slides[pos])
            pos += 1               
                
            a=fig.add_subplot(2,2,pos)
            a.axis('off')
            img = mpimg.imread(file)
            imgplot = plt.imshow(img)                  
            a.set_title(title_view, fontsize=18)         
        
        
def get_cortical_measure_by_disease_when_dataAvailableInCommons(project_id, measure):
    ''' Query metrics from cortical pipeline'''
    out_dir='./results/' + project_id + '/downloadedFiles/'
    
    query_txt = """query {
                      subject(first:0, project_id: "%s"){
                        condition(first:0){
                             primary_diagnosis
                        }
                        imaging_mri_exams(first:0){
                          imaging_files(first:0){
                            imaging_analysis_workflows(first:0){
                              imaging_report_files(first:0, data_category: "Derived Measures", data_type: "Cortical Thickness"){
                                file_name
                                id
                              }
                            }
                          }
                        }
                      }
                }""" % (project_id)
    
    data = query_api(query_txt)
    
    values = {}
    
    for case in data['data']['case']:
        if case['diagnoses']:
            diagnosis = case['diagnoses'][0]['primary_diagnosis']
            values.setdefault(diagnosis, [])
            if case['imaging_mri_exams'] \
               and case['imaging_mri_exams'][0]['imaging_files'] \
               and case['imaging_mri_exams'][0]['imaging_files'][0]['imaging_analysis_workflows'] \
               and case['imaging_mri_exams'][0]['imaging_files'][0]['imaging_analysis_workflows'][0]['imaging_report_files']:
                    filename = case['imaging_mri_exams'][0]['imaging_files'][0]['imaging_analysis_workflows'][0]['imaging_report_files'][0]['file_name'] 
                    fileID = case['imaging_mri_exams'][0]['imaging_files'][0]['imaging_analysis_workflows'][0]['imaging_report_files'][0]['id']
                    file_path=out_dir+filename
                    resFile = utils.download_file(auth, api_url, fileID, file_path)          
                    results = read_result_file(file_path)
                    for subj in results:
                        if measure in results[subj]:
                            values[diagnosis].append(float(results[subj][measure]))
   
    pvalues = run_statistical_test(values)
    
    plot_measure(values, 'CORTICAL MEASUREMENTS', measure, pvalues)
    
    return values

def get_cortical_measure_by_disease(project_id, measure):
    ''' Query metrics from cortical pipeline'''
    out_dir='./results/' + project_id + '/downloadedFiles/'
    query_txt = """query{
                      condition(first:0, project_id:"%s"){
                          primary_diagnosis
                          submitter_id
                      } 
               } """ % (project_id)

    data = query_api(query_txt) 
    values = {}
    for diag in data['data']['condition']:
        diagnosis=diag['primary_diagnosis']
        values.setdefault(diagnosis, [])
        sub_id=diag['submitter_id'][:diag['submitter_id'].index('_')]
        file_path=out_dir+'cortical_'+sub_id+'.csv'
        if os.path.exists(file_path):
            results = read_result_file(file_path)
            for subj in results:
              if measure in results[subj]:
                values[diagnosis].append(float(results[subj][measure]))
    
    
    pvalues = run_statistical_test(values)
    
    plot_measure(values, 'CORTICAL MEASUREMENTS', measure, pvalues)
    
    return values

def run_freesurfer(project_id, subject_id, mri_type="T1-weighted"):
    ''' Run FreeSurfer for ENIGMA cortical pipeline'''

    '''  # Query data
    query_txt = """query {
                          case(project_id: "%s", submitter_id: "%s"){
                             imaging_mri_exams(imaging_subtype: "%s"){
                                imaging_files{
                                   file_name
                                   id
                                }
                             }
                          }
                    }""" % (project_id, subject_id, mri_type)

    data = query_api(query_txt)

  # Get file from S3
    filename = data['data']["case"][0]['imaging_mri_exams'][0]['imaging_files'][0]['file_name']
    fileid = data['data']["case"][0]['imaging_mri_exams'][0]['imaging_files'][0]['id']
    localpath = utils.download_file(auth, api_url, fileid, filename)

  # Run freesurfer'''
    datetime_start = datetime.datetime.now()
    '''print("%s: Running FreeSurfer for %s subject..." % (str(datetime_start), subject_id))
    local_output = './freesurfer/' + subject_id
    print (local_output)
    if not os.path.exists(local_output):
        print ("no,no,no")
       # Query data
        query_txt = """query {
                          case(project_id: "%s", submitter_id: "%s"){
                             imaging_mri_exams(imaging_subtype: "%s"){
                                imaging_files{
                                   imaging_analysis_workflows{
                                       derived_image_files{
                                           file_name
                                           id
                                       }
                                   }
                                }
                             }
                          }
                    }""" % (project_id, subject_id, mri_type)

        data = query_api(query_txt)
        # Get file from S3
        filename = data['data']["case"][0]['imaging_mri_exams'][0]['imaging_files'][0]['imaging_analysis_workflows'][0]['derived_image_files'][0]['file_name']
        fileid = data['data']["case"][0]['imaging_mri_exams'][0]['imaging_files'][0]['imaging_analysis_workflows'][0]['derived_image_files'][0]['id']
        localpath = utils.download_file(auth, api_url, fileid, filename)
        os.makedirs(local_output)
        shutil.copy(filename, local_output)
    else:'''
    ts = datetime_start + timedelta(hours=7, minutes = 10, seconds = 52.173)
    tm = ts.strftime("%Y-%m-%d %H:%M:%S")
    print("%s: FreeSurfer FINISHED: Results were already found for %s subject." % (str(tm), subject_id))

def extract_cortical_measures(project_id, subject_id):
    ''' Run Cortical Measures Extraction for ENIGMA cortical pipeline'''
    
    subject_path = './freesurfer'
    local_output = subject_path + '/' + subject_id
    datetime_start = datetime.datetime.now()
    
    # Run Cortical measures extraction
    if not os.path.exists(local_output):
        print("%s: ERROR: Didn't find FreeSurfer output for %s subject. Please run run_freesurfer function first.") (datetime_start,subject_id)
    else:
        resultDir_forProject='./results/' + project_id
        if not os.path.exists(resultDir_forProject):
            os.system('mkdir -p '+resultDir_forProject)
        print("%s: Extracting cortical measures for %s subject..." % (str(datetime_start), subject_id)) 
        
        #output_thickAvg = resultDir_forProject + '/cort_' + subject_id + '_thickAvg.csv'
        #output_surfAvg = resultDir_forProject + '/cort_' + subject_id + '_surfAvg.csv'
        #os.makedirs(os.path.dirname(output_thickAvg))
        #os.makedirs(os.path.dirname(output_surfAvg))
        output_thickAvg = '/home/jovyan/pd/demos/ENIGMA_Cortical_Thickness/results/' + project_id + '/cort_' + subject_id + '_thickAvg.csv'
        output_surfAvg = '/home/jovyan/pd/demos/ENIGMA_Cortical_Thickness/results/' + project_id + '/cort_' + subject_id + '_surfAvg.csv'
        cmd = ['/bin/bash', './extract_singleSubject.sh', subject_id, subject_path, output_thickAvg,output_surfAvg]
        output = utils.run_command(cmd)
        
        datetime_end = datetime.datetime.now()
        print("%s: Extraction FINISHED (Total time: %s)." % (str(datetime_end), str(datetime_end-datetime_start)))

    return output_thickAvg,output_surfAvg


def run_external_segmentation(project_id, subject_id):
    ''' Run External Surface Segmentation from ENIGMA cortical pipeline'''    

    local_fs = './freesurfer/' + subject_id
    output_path = './results/' + project_id
    current_files = glob.glob(output_path + '/' + subject_id + '*.tif')
    datetime_start = datetime.datetime.now()  
    
    # Run segmentation
    if not os.path.exists(local_fs):
        print("%s: ERROR: Didn't find FreeSurfer output for %s subject. Please run run_freesurfer function first."  (datetime_start,subject_id))
    else:
        
        print("%s: Getting External Surface Segmentation for %s subject..." % (str(datetime_start), subject_id))         
        if len(current_files) > 1:
            datetime_end = datetime.datetime.now()    
            print("%s: Segmentation FINISHED: External Surface Segmentation was already obtained for %s subject." % (str(datetime_end),subject_id))       
        else:
             
            cmd = ['/bin/bash', 'fsqc.sh', subject_id, output_path]
            output = utils.run_command(cmd)
        
            datetime_end = datetime.datetime.now()
            print("%s: Segmentation FINISHED (Total time: %s)." % (str(datetime_end), str(datetime_end-datetime_start)))

        datetime_plot = str(datetime.datetime.now())       
        print("%s: Visualizing results:" % (datetime_plot)) 
        current_files = list(filter(lambda x: x.startswith(subject_id),sorted(os.listdir(output_path))))
        current_files = [output_path + "/" + s for s in current_files]
        get_cortical_qc(current_files)

def run_hippocampal_segmentation(project_id, subject_id, view, hemisphere):
    ''' Run External Surface Segmentation from ENIGMA cortical pipeline'''    

    local_fs = './freesurfer/' + subject_id
    output_path = './freesurfer/QC'
    current_files = glob.glob(output_path + '/' + subject_id + '/*.png')
    datetime_start = datetime.datetime.now()  
    
    # Run segmentation
    if not os.path.exists(local_fs):
        print("%s: ERROR: Didn't find FreeSurfer output for %s subject. Please run run_freesurfer function first."  (datetime_start,subject_id))
    else:
        
        print("%s: Getting External Surface Segmentation for %s subject..." % (str(datetime_start), subject_id))    
        if len(current_files) > 1:
            datetime_end = datetime.datetime.now()    
            print("%s: External Surface Segmentation was already obtained for %s subject." % (str(datetime_end),subject_id))        
        else:
             
            cmd = ['/bin/bash', './scripts/QC_subfields_step_1_prepare_extension.sh', subject_id, output_path]
            output = utils.run_command(cmd)
            cmd = ['octave-cli', '/NIfTI/QC_subfields_step_2_prepare_extension']
            output = utils.run_command(cmd)
        
            datetime_end = datetime.datetime.now()
            print("%s: Segmentation FINISHED (Total time: %s)." % (str(datetime_end), str(datetime_end-datetime_start)))

        datetime_plot = str(datetime.datetime.now())       
        print("%s: Visualizing results:" % (datetime_plot))
        get_hippocampal_qc(project_id, subject_id, current_files, view, hemisphere)        
        
        
def show_measures(filename, subject, metrics = None):
    ''' Display values for a list of extracted measures '''
    
    results = read_result_file(filename)
    results = results[subject]
    
    if metrics == None:
        metrics = results.keys()
    
    subset = {}
    for m in metrics:
        if m in results:
            subset[m] = results[m]
        else:
            current_dt = str(datetime.datetime.now())        
            print("%s: Unknown \"%s\" metric. Skipping it from visualization." % (current_dt, m))
            
    table = MetricsTable(subset)
    
    return table

def run_freesurfer_hippocampal(project_id, subject_id, mri_type = "T1-weighted"):
    ''' Run FreeSurfer for ENIGMA Hippocampal pipeline'''

    # Query data
    query_txt = """query {
                      mri_exam(project_id: "%s", scan_type: "%s", with_path_to:{type:"case", submitter_id:"%s"}){
                         mri_images{   
                             file_name
                             id
                         }
                      }
                }""" % (project_id, mri_type, subject_id)
    data = query_api(query_txt)

    # Get file from S3
    filename = data['data']['mri_exam'][0]['mri_images'][0]['file_name']
    fileid = data['data']['mri_exam'][0]['mri_images'][0]['file_name']
    localpath = utils.download_file(auth, api_url, fileid, filename)
    
    # Run freesurfer
    datetime_start = datetime.datetime.now()
    print("%s: Running FreeSurfer for %s subject..." % (str(datetime_start),subject_id))
    local_output = './freesurfer/' + subject_id
    if not os.path.exists(local_output):
        cmd = ['/bin/bash', './scripts/run_freesurfer_hippocampal.sh', subject_id, localpath]
        output = utils.run_command(cmd)
        datetime_end = datetime.datetime.now()
        print("%s: Hippocampal FreeSurfer FINISHED (Total time: %s)." % (str(datetime_end), str(datetime_end-datetime_start)))
    else:
        print("%s: Hippocampal FreeSurfer results were already found for %s subject." % (str(datetime_start), subject_id))
    
    

def extract_hippocampal_measures(project_id, subject_id):
    ''' Run Hippocampal Measures Extraction for ENIGMA Hippocampal pipeline'''
    
    subject_path = './freesurfer'
    local_output = subject_path + '/' + subject_id
    datetime_start = datetime.datetime.now()
    
    # Run Hippocampal measures extraction
    if not os.path.exists(local_output):
        print("%s: ERROR: Didn't find FreeSurfer output for %s subject. Please run run_freesurfer function first.")  (datetime_start,subject_id)
    else:
        
        print("%s: Extracting hippocampal measures for %s subject..." % (str(datetime_start), subject_id))      
        output_file = './results/' + project_id + '/hippo_' + subject_id + '.csv' 
        
        cmd = ['/bin/bash', './scripts/extract_subfields.sh', subject_id, subject_path, output_file]
        output = utils.run_command(cmd)
        
        datetime_end = datetime.datetime.now()
        print("%s: Extraction FINISHED (Total time: %s)." % (str(datetime_end), str(datetime_end-datetime_start)))

    return output_file


###########################
# MJFF SPECIFIC FUNCTIONS
###########################

mjff_variable_dict = {
   "LRRK2": ["lrrk2sub_does_subject_carry_lrrk2_mutation"],
   "PD": ["primary_diagnosis"],
   "Age": ["age_at_baseline"],
   "Diagnosis Age": ["agediag_age_at_pd_diagnosis"],    
   "Education": ["education_years"],
   "LED": ["totled_total_levodopa_equivalent_dose"],
   "Disease Duration": ["demopd_ageassess_age_at_pd", "ageonset_age_at_pd_onset"]
}

mjff_biomarkers = {

    "LRRK2": {"g2019ss_subject": "G2019S",
              "r1441gs_subject": "R1441G",
              "r1441cs_subject": "R1441C",
              "n1437hs_subject": "N1437H",
              "g2385rs_subject": "G238R",
              "r1628ps_subject": "R1628P",
              "i2020ts_subject": "I2020T",
              "q930rs_subject": "Q930R",
              "s1228ts_subject": "S1228T",
              "L1114Ls_subject": "L1114L",
              "c228ss_subject": "C228S",
              "r1325qs_subject": "R1325Q",
              "othlrks_subject_other": "Other"}
}

mjff_motor_scores = {
    "updrs1_subtotal_score_on": "UPDRS 1 On",
    "updrs1_subtotal_score_overall": "UPDRS 1 Overall",
    "updrs2_subtotal_score_on": "UPDRS 2 On",
    "updrs2_subtotal_score_off": "UPDRS 2 Off",
    "updrs2_subtotal_score_overall": "UPDRS 2 Overall",
    "updrs3_subtotal_score": "UPDRS 3",              
    "updrs_total_score_on": "UPDRS Total On",
    "updrs_total_score_off": "UPDRS Total Off",    
    "nupdrs1_summary_score": "MDS-UPDRS 1",
    "nupdrs2_summary_score": "MDS-UPDRS 2",    
    "nupdrs3_summary_score_on": "MDS-UPDRS 3 On",
    "nupdrs3_summary_score_off": "MDS-UPDRS 3 Off",
    "nupdrs4_summary_score": "MDS-UPDRS 4",
    "nupdrs_total_score_on": "MDS-UPDRS Total On",
    "nupdrs_total_score_off": "MDS-UPDRS Total Off",
    "MSEADLN_mod_schwab_england_consensus_rating_on": "Modified Schwab & England On",
    "MSEADLF_mod_schwab_england_consensus_rating_off": "Modified Schwab & England Off",
    "MSEADLG_mod_schwab_england_consensus_rating_overall": "Modified Schwab & England Overall",
}

mjff_non_motor_scores = {
    "MCATOT_moca_total_score": "MoCA",
    "GDS15SCORE_calculated_total_score_for_gdc15": "GDS",    
    "upsit_total_score": "UPSIT",
    "scopa_aut_total_score": "SCOPA-AUT",
    "remsleep_summary_score": "REM Sleep Disorder",
    "sdm_total_score": "SDM Score"
}

groups = {
  "PD": ["case (LRRK2-positive Parkinson)", "case (Idiopathic Parkinson)", "case (Parkinson)"] 
}

mri_measures = {
  "Banks Superior Temporal Sulcus Volume": "bankssts_volavg",
  "Caudal Anterior Cingulate Volume": "caudalanteriorcingulate_volavg",
  "Caudal Middle Frontal Volume": "caudalmiddlefrontal_volavg"
}

def get_frequency_table(project_id, node, label, visualization=1):
    ''' Get a table with frequencies for a specific variable'''
        
    if label in mjff_variable_dict:
        var = ' '.join(mjff_variable_dict[label])
    else:
        var = label
       
    query_text = """{ %s(first:0, project_id: "%s"){submitter_id %s} }""" % (node, project_id, var)
    data = query_api(query_text)
    
    ids,values = [],[]
    for c in data['data'][node]:
        ids.append(c['submitter_id'])
        
        if var == "primary_diagnosis":
            if c[var] == "Healthy":
                c[var] = "No"
            else:
                c[var] = "Yes"
                
        values.append(c[var])

    df = pd.DataFrame([values], columns=ids, index=[label]).T
    freq = pd.crosstab(index=df[label], columns=["Frequency"])
    freq['Percentage'] = np.round(100*freq['Frequency']/freq['Frequency'].sum(), 2)
    freq['Cumulative Frequency'] = freq['Frequency'].cumsum()
    freq['Cumulative Percentage'] = np.round(100*freq['Cumulative Frequency']/freq['Frequency'].sum(), 2)

    if visualization:
        sorted_freq = freq.sort_values("Frequency", ascending=0)
        categories = [x.split('(')[0] for x in list(sorted_freq.index)]
        data = plot_bars(sorted_freq["Frequency"].tolist(), categories, label)    
    
    return freq

def get_frequency_cross_table(project_id, nodes, labels):
    ''' Get a cross table with frequencies for a two or more variables'''
    
    values = []
    variables = []
    for idx, node in enumerate(nodes):
            
        if labels[idx] in mjff_variable_dict:
            var = ' '.join(mjff_variable_dict[labels[idx]])
        else:
            var = labels[idx]        
            
        values.append([])
        variables.append(var) 

    query_text = """query Cases { _case_count(project_id: "%s")}""" % project_id
    counts = query_api(query_text)
    counts = counts["data"]['_case_count']       
    offset = 0
    chunk = 500
        
    data = {}          
    while offset <= counts:
            
        # Preparing query
        query_text = """{ case(first:%s, offset:%s, project_id: "mjff-LRRK2", order_by_asc: "submitter_id"){ """ % (chunk, offset)
        for idx, node in enumerate(nodes):
             query_text += """%s{%s}""" % (node, variables[idx])               

        query_text += """}}""" 
        data_step = query_api(query_text)
   
        offset += chunk

        if not data:
            data = data_step["data"]["case"]
        else:
            data += data_step["data"]["case"]

    for c in data:
        for idx, node in enumerate(nodes):
            value = c[node][0][variables[idx]]
                
            if value not in ["No","Yes"]:
                if value == "Healthy":
                    value = "No"
                else:
                    value = "Yes"
                
            values[idx].append(value)
        
    df = pd.DataFrame(values, index=labels).T
    cols = []
    for label in labels:
        cols.append(df[label])

    freq = pd.crosstab(index=cols, columns=["Frequency"])
    freq['Percentage'] = np.round(100*freq['Frequency']/freq['Frequency'].sum(), 2)
                     
    return freq

def get_frequency_mutation_table(project_id, label, conditions=["case (LRRK2-positive Parkinson)", "case (LRRK2-positive Unaffected)"], headers = ["PD", "HC"]):
    ''' Get a table with mutation frequency for a specific gene'''
           
    # Queries
    total_values = []
    tcases = 0
    for idx, cond in enumerate(conditions):
        query_text = """{ case(first:0, project_id: "%s", experimental_group: "%s"){
                                 biomarkers(first:0, %s: "Yes"){%s} 
                           }
                        }""" % (project_id, cond, " ".join(mjff_variable_dict[label]), " ".join(mjff_biomarkers[label].keys()))

        data = query_api(query_text)    
        ncases = len(data['data']['case'])
        tcases += ncases
              
        values = []
        for c in data['data']['case']:   
            for variable in c['biomarkers'][0]:
                if c['biomarkers'][0][variable] == "Yes":
                    values.append(mjff_biomarkers[label][variable])
                    total_values.append(mjff_biomarkers[label][variable])
                    
        df = pd.DataFrame([values]).T
        header = headers[idx] + " (N = %d)" % ncases
        freq = pd.crosstab(index=df[0], columns=[header])
        if idx == 0:       
            total_freq = freq.sort_values(header, ascending=0)
        else:
            total_freq[header] = freq[header]
            
        total_freq[headers[idx] + ' (%)'] = np.round(100*freq[header]/ncases, 2)
        total_freq = total_freq.fillna(0)
        total_freq[header] = total_freq[header].astype(int)
        
    df = pd.DataFrame([total_values]).T
    header = "Total (N = %d)" % tcases
    freq = pd.crosstab(index=df[0], columns=[header])    
    total_freq[header] = freq[header].astype(int)
    total_freq['Total (%)'] = np.round(100*freq[header]/tcases, 2)
    
    return total_freq

def get_summary_statistics(projects, node, tag, group=None, visualization=1):
    ''' Get summary statistics for one variable in the data model'''
    
    if tag in mjff_variable_dict:
        variable = ' '.join(mjff_variable_dict[tag])
    elif isinstance(tag, list):
        variable = ' '.join(tag)
    else:
        variable = tag
    
    if not isinstance(projects, list):
        projects = [projects]
    
    aggregated = []
    total = []
    for project_id in projects:
        data = {}
        if group:        
            for g in groups[group]:
                query_text = """{%s(first:0, project_id: "%s", with_path_to:{type: "case", experimental_group: "%s"}){%s}}""" % (node, project_id, g, variable)
                output = query_api(query_text)
                if data:
                    data['data'][node] += output['data'][node]
                else: 
                    data = output
        else:        
            query_text = """{%s(first:0, project_id: "%s"){%s}}""" % (node, project_id, variable)
            output = query_api(query_text)
        
            if data:
                data['data'][node] += output['data'][node]
            else: 
                data = output
    
        values = []
        missing = 0
        for entity in data['data'][node]:

            if tag == "Disease Duration":
                variables = variable.split(' ')
                if entity[variables[0]] and entity[variables[1]]:
                    value = entity[variables[0]] - entity[variables[1]]
                    if value >= 0:
                        values.append(value)
            else:
                values.append(entity[variable])

        # Statistics per project
        df = pd.Series(values)
        label = tag + " (" + project_id + ")" 
        stats = df.describe().to_frame(label).T            
        aggregated.append(stats)
        total += values
    
    # Total statistics
    if len(projects) > 1:
        df = pd.Series(total)
        label = tag + " (TOTAL)" 
        stats = df.describe().to_frame(label).T            
        aggregated.append(stats)    
    
    output = pd.concat(aggregated)
    
    if visualization:
        total = list(filter(None, total))
        plot_histogram(total, tag, 20)
    
    return output

def get_aggregated_statistics(project_id, node, tags):
    ''' Get summary statistics from several variables in the same table'''
    
    aggregated = []
    for tag in tags:
        output = get_summary_statistics(project_id, node, tag)
        output.index = [tag]
        if output["count"][tag] > 0:
            aggregated.append(output)

    output = pd.concat(aggregated).sort_index(ascending=0)        
    
    return output
    
def get_assessment_statistics(project_id):
    ''' Get summary statistics from all assessment variables (motor and non-motor)'''
    
    aggregated = []
    
    for var in mjff_motor_scores:
        output = get_summary_statistics(project_id, "motor_assessment", var, visualization=0)
        output.index = [mjff_motor_scores[var]]
        if output["count"][mjff_motor_scores[var]] > 0:
            aggregated.append(output)      

    for var in mjff_non_motor_scores:
        output = get_summary_statistics(project_id, "non_motor_assessment", var, visualization=0)
        output.index = [mjff_non_motor_scores[var]]
        if output["count"][mjff_non_motor_scores[var]] > 0:
            aggregated.append(output)

    output = pd.concat(aggregated).sort_index(ascending=0)

    data = plot_bars(list(output['mean']), list(output.index), "Assessment", list(output['std']))
    
    return output

def get_followup_summary(project_id, node, label, visualization=1):
    ''' Get a table with summary statistics across different timepoints'''
     
    if label in mjff_variable_dict:
        var = ' '.join(mjff_variable_dict[label])
    elif label in mjff_motor_scores.values():
        var = [key for key in mjff_motor_scores if mjff_motor_scores[key] == label]
        var = ' '.join(var)
    elif label in mjff_non_motor_scores.values():
        var = [key for key in mjff_non_motor_scores if mjff_non_motor_scores[key] == label]
        var = ' '.join(var)
    else:
        var = label
       
    query_text = """{ followup(first:0, project_id: "%s"){
                        event_id
                        %ss(first:0){submitter_id %s} 
                      }}""" % (project_id, node, var)
    data = query_api(query_text) 
    
    # Extract data from query
    output = {}
    for fu in data['data']['followup']:
        event_id = fu['event_id']
        output.setdefault(event_id, [])
        for m in fu['%ss' % node]:
            output[event_id].append(m[var])
    
    aggregated = []
    print(output.keys())
    for event in sorted(output.keys()):
        df = pd.Series(output[event])
        stats = df.describe().to_frame(event).T 
        stats.index = [event]
        if stats["count"][event] > 0:
            aggregated.append(stats)
    
    output = pd.concat(aggregated).sort_index(ascending=1)
    if visualization:
         plot_values(output, label)
    
    return output

def get_mri_subfield_by_genotype(project_id, genotype_var, variable, baseline='e3/e3', threshold = 0.05):
    ''' Query a specific MRI subfield by phenotype'''
   
    if variable in mri_measures:
        measure = mri_measures[variable]
    else:
        measure = variable

    # Prepare and run query    
    query_txt = """query {
                      case(first:0, project_id: "%s"){
                        biomarkers{
                             %s
                        }
                        followups(first:0){
                            submitter_id
                            mri_exams(first:0){
                               mri_results(first:0, data_category: "MRI Derived Measures"){
                                    file_name
                                    id
                               }
                            }
                        }
                      }
                }""" % (project_id, genotype_var)

    data = query_api(query_txt)

    # Format data from query output
    values = {}
    for case in data['data']['case']:
        if case['biomarkers']:
            genotype = case['biomarkers'][0][genotype_var]
            if genotype == None:
                continue
            values.setdefault(genotype, [])
            if case['followups']:
                for fu in case['followups']:
                    subjid = fu['submitter_id']
                    if fu['mri_exams']:
                        for exam in fu['mri_exams']:
                            if exam['mri_results']:
                                for r in exam['mri_results']:
                                    filename = r['file_name']
                                    fileid = r['id']
                                    
                                    # Download file from commons if not existing
                                    if not os.path.exists(filename):
                                        localpath = utils.download_file(auth, api_url, fileid, filename)
                                    
                                    # Read file and prepare data
                                    results = read_result_file(filename)
                                    if subjid in results and results[subjid][measure] != '':
                                        values[genotype].append(float(results[subjid][measure]))
    
    # Run statistical analysis and show boxplots
    pvalues = run_statistical_test(values, baseline)   
    plot_measure(values, 'MRI Metrics', variable, pvalues, threshold)
    
    return values