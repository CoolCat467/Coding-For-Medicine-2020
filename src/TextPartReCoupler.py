"""Couples overlapping sequences of text (and DNA) back together after sequencing."""

# Programmed by CoolCat467

from __future__ import annotations

__title__ = "Recoupler"
__version__ = "1.0.0"


def recouple(lst: list[str], scan_length: int = 8, start: str = "") -> str:
    """Program that will re-couple text from a list that has been broken into many overlapping parts."""
    # Get all the separate chunks of data in a list, no repeats.
    data = sorted(set(lst))
    if scan_length < 4:
        scan_length = min(round(min([len(i) for i in data]) / 2), 4)
        print(f"{scan_length = }")
    # Get all the words into a dictionary
    words_list = {}
    beginnings = {}
    # For each thing to test in the dictionary
    for test in data:
        # If we haven't seen this word before, set it as a possible
        # beginning value
        if test[0:scan_length] not in beginnings:
            beginnings[test[0:scan_length]] = True
        # For each index point our test string,
        for i in range(len(test) - scan_length):
            # Get the word it could be
            word = test[i : i + scan_length]
            # Add the word to the words dict, pointing to the last letter of
            # our test string.
            words_list[word] = test[i + scan_length]
            # If it's not the beginning
            if i > 0:
                beginnings[word] = False
    # If the starting word was not defined,
    if not start:
        # Set it to a random valid word
        for word in beginnings:
            if beginnings[word]:
                test = str(word)
                break
    else:
        # Otherwise, try to use it
        test = str(start)
        yes = False
        # If it's invalid, break
        if test not in words_list:
            for word in words_list:
                if test in word:
                    test = word
                    yes = True
                    break
        if not yes:
            raise ValueError(
                "Start word argument is invalid, does not exist in word list with scan_length. Try having start word len be equal to scan_length argument number.",
            )
    # Free up some memory
    del beginnings
    # Get data that points to each other
    string = test
    count = len(words_list.keys()) * scan_length
    while count:
        # If the test word is valid
        if test in words_list:
            # Add it to the string
            string += words_list[test]
            test = test[1:scan_length] + words_list[test]
            count -= 1
        else:  # Otherwise, decrement scan by one
            count = 0
    ##    # Once we have data that points to other bits, put them together properly
    ##    data = str(string[0])
    ##    # For each bit of data to add,
    ##    for i in string[1:]:
    ##        # For all possible index positions in the new bit of data,
    ##        for ii in range(len(i)):
    ##            # Get the current data's tail with ii
    ##            piece = data[-ii:]
    ##            # If the tail of this index is the start of this bit of data,
    ##            if piece == i[:len(piece)]:
    ##                # Combine the two chucks properly and stop checking indexes for this new chunk.
    ##                data = data[:len(data)-ii] + i
    ##                break
    return string


def run():
    """Recouple data."""
    print("Running...")
    data = [
        "ossible-copyin",
        "tely-suggests-a",
        "our-notice-that",
        "ice-that-the-sp",
        "ible-copying-me",
        "ggests-a-possi",
        "-suggests-a-pos",
        "anism-for-the-g",
        "lated-immediate",
        "ped-our-notice-",
        "y-suggests-a-p",
        "e-copying-mecha",
        "ying-mechanism-",
        "ts-a-possible-c",
        "ossible-copying",
        "ying-mechanism",
        "the-genetic-mat",
        "-immediately-su",
        "pecific-pairing",
        "uggests-a-possi",
        "-copying-mecha",
        "copying-mechani",
        "aped-our-notice",
        "immediately-sug",
        "ted-immediatel",
        "ng-we-have-post",
        "ing-mechanism-f",
        "mediately-sugge",
        "iately-suggests",
        "uggests-a-poss",
        "ur-notice-that-",
        "-have-postulate",
        "fic-pairing-we-",
        "we-have-postul",
        "e-genetic-mater",
        "ostulated-immed",
        "stulated-immed",
        "ed-immediately-",
        "ce-that-the-spe",
        "notice-that-the",
        "-pairing-we-hav",
        "ng-we-have-pos",
        "-for-the-geneti",
        "as-not-escaped-",
        "t-has-not-escap",
        "-postulated-imm",
        "-postulated-im",
        "t-the-specific-",
        "ests-a-possible",
        "-not-escaped-ou",
        "-that-the-spec",
        "otice-that-the-",
        "pairing-we-hav",
        "has-not-escaped",
        "e-have-postulat",
        "ely-suggests-a-",
        "-specific-pair",
        "ble-copying-mec",
        "-genetic-materi",
        "-our-notice-tha",
        "mmediately-sugg",
        "caped-our-notic",
        "tice-that-the-s",
        "that-the-speci",
        "-copying-mechan",
        "possible-copyin",
        "ve-postulated-i",
        "ed-immediately",
        "ng-mechanism-fo",
        "specific-pairin",
        "e-specific-pair",
        "ce-that-the-sp",
        "ulated-immediat",
        "-the-specific-p",
        "escaped-our-not",
        "diately-sugges",
        "at-the-specific",
        "-that-the-speci",
        "gests-a-possib",
        "ecific-pairing-",
        "not-escaped-our",
        "sm-for-the-gene",
        "or-the-genetic",
        "ur-notice-that",
        "ing-we-have-pos",
        "ave-postulated-",
        "stulated-immedi",
        "aped-our-notic",
        "e-that-the-spec",
        "a-possible-copy",
        "he-specific-pai",
        "-have-postulat",
        "ic-pairing-we-h",
        "-possible-copyi",
        "It-has-not-esca",
        "-the-genetic-m",
        "g-mechanism-for",
        "sible-copying-m",
        "hat-the-specif",
        "r-notice-that-",
        "-notice-that-th",
        "e-postulated-i",
        "d-our-notice-th",
        "cific-pairing-w",
        "iately-suggest",
        "pairing-we-have",
        "tulated-immedi",
        "e-postulated-im",
        "r-the-genetic-",
        "we-have-postula",
        "r-notice-that-t",
        "-escaped-our-n",
        "chanism-for-the",
        "t-escaped-our-n",
        "pecific-pairin",
        "hanism-for-the-",
        "for-the-genetic",
        "ly-suggests-a-p",
        "-has-not-escap",
        "r-the-genetic-m",
        "he-genetic-mat",
        "c-pairing-we-ha",
        "-the-genetic-ma",
        "diately-suggest",
        "pying-mechanism",
        "for-the-geneti",
        "g-we-have-postu",
        "y-suggests-a-po",
        "at-the-specifi",
        "scaped-our-noti",
        "s-not-escaped-",
        "tulated-immedia",
        "e-that-the-spe",
        "ism-for-the-gen",
        "nism-for-the-g",
        "-for-the-genet",
        "postulated-imme",
        "ately-suggests-",
        "ssible-copying-",
        "opying-mechanis",
        "c-material",
    ]
    print('"', recouple(data), '"')
    ##    data = {len(i):i for i in [recouple(DATA, DATA[i][:6], 6) for i in range(len(DATA))]}
    ##    print(max(data.keys()), data[max(data.keys())])
    print("Done")


if __name__ == "__main__":
    run()
