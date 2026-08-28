#!/usr/bin/env python3
resolved={("BRXN_ALK_B","GBFBEFJJ_01512")}
hits=[
 ("BRXN_ALK_B","GBFBEFJJ_01512"),
 ("BRXN_OTHER","GBFBEFJJ_01512"),
 ("BRXN_ALK_B","OTHER_LOCUS"),
]
kept=[x for x in hits if x not in resolved]
assert ("BRXN_ALK_B","GBFBEFJJ_01512") not in kept
assert ("BRXN_OTHER","GBFBEFJJ_01512") in kept
assert ("BRXN_ALK_B","OTHER_LOCUS") in kept
print("PASS pair-specific skip behavior")
