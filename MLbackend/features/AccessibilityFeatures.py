import tempfile
from pathlib import Path
from subprocess import call
import os
from features.Features import Features
import subprocess
import requests
import re




class AccessibilityFeatures(Features):
    def extract_features(self):
        self._features_dict["accessibility"] = self.accessibility()

    def accessibility(self):  # 37*10 = 370
        url = "http://localhost:3031/wsl_api/plfold?region_seq=" + self._region_sequence
        response = requests.get(url)
        ACC_allstr = response.text
        ACC_allstr = ACC_allstr.split("\\n")
        ACC_allstr = ACC_allstr[: len(ACC_allstr) - 1]
        ACC_allstr = [s + '\\n' for s in ACC_allstr] 
        ACC_allstr = [s.lstrip('","') for s in ACC_allstr]
        ACC_allstr[0] = ACC_allstr[0][2:]
        acc_score_matrix = []
        for line in ACC_allstr[2:]:
            l = line.strip().replace('NA', str(0.97)).split()
            acc_score_matrix.append(l)
        acc_score_matrix = [[0.0] * 11] * 37 + acc_score_matrix + [[0.0] * 11] * 37
        acc_score_matrix = [line[0].split("\\t") if isinstance(line[0], str) else line for line in acc_score_matrix]
        acc_score_matrix_segment = acc_score_matrix[self._end + 15:self._end + 52]
        acc_score_matrix_segment = [line.split("\\t") if isinstance(line, str) else line for line in acc_score_matrix_segment]
        for i in range(len(acc_score_matrix_segment)):
            for j in range(len(acc_score_matrix_segment[i])):
                if '\\n' in acc_score_matrix_segment[i][j]:
                    acc_score_matrix_segment[i][j] = acc_score_matrix_segment[i][j].replace('\\n', "")
        ACC = {}
        try:
            for i in range(1, 38):
                for j in range(1, 11):
                    key = 'Acc_P%s_%sth' % (str(i), str(j))
                    ACC[key] = float(acc_score_matrix_segment[i - 1][j])
        except:
            print('was except')
            print(self._region_sequence)
            ACC[key] = float(0.0)
        return ACC
