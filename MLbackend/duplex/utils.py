def mix_inter_bulge_char(a, b):
    assert a==' ' or b==' ', f"mix_inter_bulge error: both interaction and bulge have nt value. " \
                             f"one of them should be empty (space)." \
                             f"a={a}   b={b}"
    return chr(ord(a) + ord(b) - ord(' '))


def mix_inter_bulge_seq(s1, s2):
    r=""
    for i in range(max(len(s1), len(s2))):

        try:
            a = s1[i]
        except IndexError:
            a = " "

        try:
            b = s2[i]
        except IndexError:
            b = " "

        r += mix_inter_bulge_char(a, b)
    return r.replace(" ", "")

