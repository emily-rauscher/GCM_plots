;REMEMBER TO COPY APPROPRIATE FORT.2 AND FORT.26 FILES

;grab needed values from fort.2 and save to .txt files
;$grabber
$grabberN18

nl = 18
oom = 1.2
unstable = 1

teq = rd_tfile('teq.txt', nl, /convert)
;zheight = rd_tfile('zheight.txt', nl, /convert)
dtep = rd_tfile('dtep.txt', 1, /convert)
tgr = rd_tfile('tgr.txt', 1, /convert)
ztrop = rd_tfile('ztrop.txt', 1, /convert)
alr = rd_tfile('alr.txt', 1, /convert)
ga = rd_tfile('ga.txt', 1, /convert)
gascon = rd_tfile('gascon.txt', 1, /convert)

sigtrop = ((tgr-ztrop*alr)/tgr)^(ga/alr/gascon)

sigtrop = sigtrop[0]
dtep = dtep[0]

Tprofiles, oom = oom, teq = teq, sigtrop = sigtrop, dtep = dtep,  xmin = xmin, xmax = xmax, $
  zheight = zheight, oldbug = oldbug, unstable = unstable

delvarx, xmin, xmax

