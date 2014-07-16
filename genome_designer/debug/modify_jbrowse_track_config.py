"""Functions for manipulating JBrowse configs.

NOTE: User responsible for managing backups / not breaking anything.
"""

import json

TRACK_LIST_CONFIG = '/dep_data/temp_data/projects/3bc32fc9/ref_genomes/01166f51/jbrowse/trackList.json'

def main():
  with open(TRACK_LIST_CONFIG) as fh:
    config_json = json.loads(fh.read())

  tracks = config_json['tracks']
  for track in tracks:
    track['chunkSizeLimit'] = 1000000000
    track['maxHeight'] = 10000


  with open(TRACK_LIST_CONFIG, 'w') as output_fh:
    output_fh.write(json.dumps(config_json))

if __name__ == '__main__':
  main()
