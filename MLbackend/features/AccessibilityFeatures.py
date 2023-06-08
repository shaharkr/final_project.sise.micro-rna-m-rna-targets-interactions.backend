import tempfile
from pathlib import Path
from subprocess import call

from features.Features import Features


class AccessibilityFeatures(Features):
    def extract_features(self):
        self._features_dict["accessibility"] = self.accessibility()

    def accessibility(self):  # 37*10 = 370
        with tempfile.TemporaryDirectory() as tmpdirname:
            acc_file: Path = Path(tmpdirname) / "mrna_acc.fa"
            with open(acc_file, 'w') as f:
                f.write(self._region_sequence.replace('-', '') + '\n')

            cmd = 'RNAplfold -W 80 -L 40 -u 10 < {}'.format(acc_file)
            status = call(cmd, cwd=Path(tmpdirname).resolve(), shell=True)

            acc_outfile: Path = Path(tmpdirname) / "plfold_lunp"
            with acc_outfile.open('r') as f:
                ACC_allstr = f.readlines()

        acc_score_matrix = []
        for line in ACC_allstr[2:]:
            l = line.strip().replace('NA', str(0.97)).split()
            acc_score_matrix.append(l)

        acc_score_matrix = [[0.0] * 11] * 37 + acc_score_matrix + [[0.0] * 11] * 37
        acc_score_matrix_segment = acc_score_matrix[self._end + 15:self._end + 52]
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
