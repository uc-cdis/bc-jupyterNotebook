import json
import requests
import shutil
    
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
