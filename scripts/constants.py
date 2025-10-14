from dataclasses import dataclass
from enum import Enum

@dataclass
class Language:
    id: str # ISO alpha-2 code
    name: str
    sub_family: str
    family: str
    


LANGUAGES = [
    # Indo-European
    Language("en", "English", "Germanic", "Indo-European"),
    Language("de", "German", "Germanic", "Indo-European"),
    Language("sv", "Swedish", "Germanic", "Indo-European"),
    Language("is", "Icelandic", "Germanic", "Indo-European"),
    Language("it", "Italian", "Romance", "Indo-European"),
    Language("fr", "French", "Romance", "Indo-European"),
    Language("es", "Spanish", "Romance", "Indo-European"),
    Language("gl", "Galician", "Romance", "Indo-European"),
    Language("pt", "Portuguese", "Romance", "Indo-European"),
    Language("cs", "Czech", "Slavic", "Indo-European"),
    Language("pl", "Polish", "Slavic", "Indo-European"),
    Language("ru", "Russian", "Slavic", "Indo-European"),
    Language("hi", "Hindi", "Indo-Aryan", "Indo-European"),

    # Other families
    Language("ar", "Arabic", "Semitic", "Afro-Asiatic"),
    Language("zh", "Chinese", "Sinitic", "Sino-Tibetan"),
    Language("fi", "Finnish", "Finnic", "Uralic"),
    Language("id", "Indonesian", "Malayo-Polynesian", "Austronesian"),
    Language("ja", "Japanese", "Japonic", "Japonic"),
    Language("ko", "Korean", "Koreanic", "Koreanic"),
    Language("th", "Thai", "Tai", "Kra-Dai"),
    Language("tr", "Turkish", "Oghuz", "Turkic"),
]


class ConlluColumns:
    ID = 0
    FORM = 1
    LEMMA = 2
    UPOS = 3
    XPOS = 4
    FEATS = 5
    HEAD = 6
    DEPREL = 7
    DEPS = 8
    MISC = 9
    TOTAL = 10


pass
