import subprocess
import glob
import os
import json
import requests
import shutil

def run_command(cmd):
    ''' Run any command for pipeline'''
    
    try:
        output = subprocess.check_output(cmd, stderr=subprocess.STDOUT)            
    except Exception as e:
        output = e.output
        output = output.decode('utf-8')
        print ("ERROR:" + output)

    return output

def get_file_from_s3(fname, project_id):
    ''' Download file from S3'''  
   
    # Create folder
    local_path = './results/' + project_id
    file_path = local_path + '/' + fname
    if not os.path.exists(local_path):
        os.makedirs(local_path)
    if os.path.exists(file_path):
        return file_path
    
    # Get bucket name and path
    s3_path = 's3://bhc-data/' + project_id + '/' + fname

    proxy = 'http://cloud-proxy.internal.io:3128'
    os.environ['http_proxy'] = proxy
    os.environ['https_proxy'] = proxy 

    # Getting files        
    cmd = ['aws', 's3', 'cp', s3_path, file_path, '--profile', 'bhc']
    exit = run_command(cmd)
    
    del os.environ['http_proxy']
    del os.environ['https_proxy']       
        
    return file_path

def download_file(auth, api_url, fileid, output):

    download_url = api_url + 'user/data/download/'
    file_url = requests.get(download_url+fileid, headers={'Authorization': 'bearer '+ auth.json()['access_token']}).text
    file_url = json.loads(file_url)

    if 'url' in file_url:
        response = requests.get(file_url['url'], stream=True)
        with open(output, 'wb') as out_file:
            shutil.copyfileobj(response.raw, out_file)
    else:
        print(file_url)
        response = file_url

    return response
