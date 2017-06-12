"""
Pretty much does what e-fetch does

www.youtube.com/watch?v=E7wB__M9fdw 
(Intro to python web scraping by Chris Reeves)

"""
root = 'http://www.ncbi.nlm.nih.gov'

def build_template(page,term):
    template = open(page,'wb')
    template.write(b'<!-- Latest compiled and minified CSS -->')
    template.write(b'<link rel="stylesheet" href="http://maxcdn.bootstrapcdn.com/bootstrap/3.2.0/css/bootstrap.min.css">')
    template.write(b'<!-- Optional theme -->')
    template.write(b'<link rel="stylesheet" href="http://maxcdn.bootstrapcdn.com/bootstrap/3.2.0/css/bootstrap-theme.min.css">')
    template.write(b'<!-- Latest compiled and minified JavaScript -->')
    template.write(b'<script src="http://maxcdn.bootstrapcdn.com/bootstrap/3.2.0/js/bootstrap.min.js"></script>')
    for title in list_title(term):
        template.write(bytes(title,'UTF-8'))
        template.write(b"<br>")
        template.write(bytes(get_abstract(find_link(title)[0]),'UTF-8'))
        
def list_title(term,with_link=True):
    from urllib.request import urlopen
    import re
    
    titles = []
    htmlfile = urlopen('http://www.ncbi.nlm.nih.gov/pubmed/?term={0}'.format(term))
    encoding = htmlfile.headers.get_content_charset()

    htmltext = htmlfile.read().decode(encoding)

    regex = '<p class="title"><.+?></p>'
    link_regex = '<a\s+(?:[^>]*?\s+)?href="([^"]*)"'
    title_regex = 'linksrc=docsum_title">.+?</a>'
    
    pattern = re.compile(regex)
    link_pattern = re.compile(link_regex)
    title_pattern = re.compile(title_regex)
    
    pubs = re.findall(pattern,htmltext)
    for pub in pubs:
        lmatch = re.search(link_pattern, pub)
        tmatch = re.search(title_pattern, pub)
        title = (pub[tmatch.start()+22:tmatch.end()-4]).encode(encoding)
        if(with_link == True): 
            link = full_link(root,pub[lmatch.start()+9:lmatch.end()-1])
            titles.append('<a href={0}>{1}</a>'.format(link,title.decode(encoding)))
        else:
            titles.append(title.decode(encoding))
    return(titles)

# usage: print(get_abstract("http://www.ncbi.nlm.nih.gov/pubmed/25090088",1))
def get_abstract(link,style=True):
    from urllib.request import urlopen
    import re
    
    htmlfile = urlopen('http://'+link)
    encoding = htmlfile.headers.get_content_charset()
    htmltext = htmlfile.read().decode(encoding)

    regex = '<div class="abstr"><.+?></div>'
    iso_regex = '<\s*p[^>]*>([^<]*)<\s*\/\s*p\s*>'
    
    pattern = re.compile(regex)
    iso_pattern = re.compile(iso_regex)
    try:
        abstract_all = re.findall(pattern,htmltext)[0]
        abstract = re.split(iso_pattern,abstract_all)[1]
    except IndexError:
        abstract_all = abstract = ''
    return abstract_all if style is True else abstract

def full_link(root,rel):
    return root + rel

def find_link(line):
    import re
    url_regex = '(?:[*]*?\s+)?http://([^>]*)'
    url_pattern = re.compile(url_regex)
    url = re.findall(url_pattern,line) # should only be one though...
    return url

from Bio import Entrez, Medline
from datetime import datetime

# Make sure you change this to your email
Entrez.email = 'brianfbb@yahoo.com'

def fetch(t, s):
    h = Entrez.esearch(db='pubmed', term=t, retmax=5, sort=s)
    idList = Entrez.read(h)['IdList']

    if idList:
        handle = Entrez.efetch(db='pubmed', id=idList,
                               rettype='medline', retmode='text')
        records = Medline.parse(handle)

        for record in records:
            title = record['TI']
            author = ', '.join(record['AU'])
            source = record['SO']
            pub_date = datetime.strptime(record['DA'], '%Y%m%d').date()
            pmid = record['PMID']

            print("Title: %s\nAuthor(s): %s\nSource: %s\n"\
                    "Publication Date: %s\nPMID: %s\n" % (title, author,
                        source, pub_date, pmid))

print('-- Sort by publication date --\n')
fetch('Dmel wings', 'pub date')

print('-- Sort by first author --\n')
fetch('Dmel wings', 'author')

