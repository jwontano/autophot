'''
John W. Montano
LCO API program to see which proposals I am a part of and I have options to download images to my station. A lot of this code is based on the LCO API documentation and github. 
https://github.com/LCOGT/
'''
import argparse, requests, os, time, sys
from datetime import datetime, timedelta
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-startdate', default=None, help='Start date for search range [YYYY-MM-DD]')
    parser.add_argument('-enddate', default=None, help='End date for search range [YYYY-MM-DD]')
    parser.add_argument('-proposal_id', default=None, help='Proposal id to search with')
    parser.add_argument('-target', default=None, help='Object target')
    parser.add_argument('-list', default=False, help='Will list proposals')
    parser.add_argument('-public', default=False, help='Search public proposals if True')
    args = parser.parse_args()
    return args


def get_token(plist=False):
    response = requests.post(
        'https://archive-api.lco.global/api-token-auth/',
        data = {
            'username': 'korbinite5646',
            'password': 'd1ckbutt'
        }
    ).json()
    token = response.get('token')
    headers = {'Authorization': 'Token ' + token}
    if plist:
        profile_response = requests.get('https://archive-api.lco.global/profile/',headers=headers)
        try: 
            profile_response.raise_for_status()
        except requests.exceptions.HTTPError as exc:
            print('Request failed: {}'.format(profile_response.content))
            raise exc
        profile_dict = profile_response.json()
        print('Username: {} proposals'.format(profile_dict['username']))
        print('You are a member of {} proposals'.format(len(profile_dict['profile']['proposals'])))
        for p in profile_dict['profile']['proposals']:
            print(p)
    return token


def download(start_date,end_date,proposal_id,data_directory,target=None,public=False):

    token = get_token()
    headers = {'Authorization': 'Token ' + token}
    if target == None:
        archive_response = requests.get('https://archive-api.lco.global/frames/?reduction_level=91&limit=50&proposal_id='+proposal_id+'&'+'start='+start_date+'&'+'end='+end_date+' 23%3A59&&',headers=headers).json()
    elif public:
        archive_response = requests.get('https://archive-api.lco.global/frames/?reduction_level=91&limit=100&'+'proposal_id='+proposal_id+'&'+'start='+start_date+'&'+'end='+end_date+' 23%3A59&'+'target_name='+target+'&public=true',headers=headers).json()
    else:
        archive_response = requests.get('https://archive-api.lco.global/frames/?reduction_level=91&limit=20&'+'proposal_id='+proposal_id+'&'+'start='+start_date+'&'+'end='+end_date+' 23%3A59&'+'target_name='+target+'&public=false',headers=headers).json()
    frames = archive_response['results']
    files = []
    for i,j,k in os.walk('/home/korbinite5646/AGN_home/FITS/'):
        if k != []:
            files = files + k
        else:
            continue
    while True:
        for frame in frames:
            if (str(frame['filename'][:-3]) in files) or (str(frame['filename']) in files):
                print('{} file exists already skipping...'.format(frame['filename']))
                continue
            else:
                with open(os.path.join(data_directory, frame['filename']), 'wb') as f:
                    print('Downloading {}...\n'.format(frame['filename']))
                    f.write(requests.get(frame['url']).content)
        if archive_response.get('next'):
            archive_response = requests.get(archive_response['next'], headers=headers).json()
            frames = archive_response['results']
        else:
            break

if __name__ == '__main__':
    args = parse_args()

    start_date = args.startdate
    end_date = args.enddate
    propID = args.proposal_id
    target = args.target
    prop_list = args.list
    public = args.public

    if propID == None and prop_list == False:
        print('Please input a proposal identification to choose from. If you would like to see the proposal list use the -list option.')
        exit()
    if prop_list:
        token = get_token(True)
        exit()
    data_d = '/home/korbinite5646/AGN_home/FITS/API/'
    if end_date is None and start_date is None:
        end_date =  datetime.today()
        start_date = end_date - timedelta(days=7)
        end_date = end_date.strftime("%Y-%m-%d")
        start_date = start_date.strftime("%Y-%m-%d")
        print('Searching last 7 days from {} to {}...\n'.format(start_date,end_date))
    elif end_date is None:
        end_date = datetime.today()
        end_date = end_date.strftime("%Y-%m-%d")
        print('Searching date range from {} to {}...\n'.format(start_date,end_date))
        
    download(start_date,end_date,propID,data_d,target,public=public)
