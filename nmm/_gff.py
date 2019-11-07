from typing import NamedTuple, List, Union

Item = NamedTuple(
    "Item",
    [
        ("seqid", str),
        ("source", str),
        ("type", str),
        ("start", int),
        ("end", int),
        ("score", Union[float, str]),
        ("strand", str),
        ("phase", Union[int, str]),
        ("attributes", str),
    ],
)


class GFFWriter:
    def __init__(self):
        self._items: List[Item] = []

    def append(self, item: Item):
        self._items.append(item)

    def dump(self, fp):
        fp.write("##gff-version 3\n")
        for item in self._items:
            cols = [
                item.seqid,
                item.source,
                item.type,
                str(item.start),
                str(item.end),
                str(item.score),
                item.strand,
                str(item.phase),
                item.attributes,
            ]
            fp.write("\t".join(cols))
            fp.write("\n")


# gff-version 3
# ctg123  .  exon  1300  1500  .  +  .  ID=exon00001
# ctg123  .  exon  1050  1500  .  +  .  ID=exon00002
# ctg123  .  exon  3000  3902  .  +  .  ID=exon00003
# ctg123  .  exon  5000  5500  .  +  .  ID=exon00004
# ctg123  .  exon  7000  9000  .  +  .  ID=exon00005
