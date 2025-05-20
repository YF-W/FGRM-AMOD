# FGRM-AMOD
FGRM-AMOD: An Adaptive Multi-View Outlier Detection Algorithm based on Fuzzy Rough Set Multi-Granularity(2024,Code)

<table class="MsoTable15Plain3" border="1" cellspacing="0" cellpadding="0" align="left" width="589" style="width:442.0pt;border-collapse:collapse;border:none;
 mso-border-alt:solid windowtext .5pt;mso-table-overlap:never;mso-yfti-tbllook:
 1184;mso-table-lspace:9.0pt;margin-left:6.75pt;mso-table-rspace:9.0pt;
 margin-right:6.75pt;mso-table-anchor-vertical:paragraph;mso-table-anchor-horizontal:
 margin;mso-table-left:left;mso-table-top:.05pt;mso-padding-alt:0cm 5.4pt 0cm 5.4pt;
 mso-border-insideh:.5pt solid windowtext;mso-border-insidev:.5pt solid windowtext">
 <tbody><tr style="mso-yfti-irow:-1;mso-yfti-firstrow:yes;mso-yfti-lastfirstrow:yes;
  height:16.15pt">
  <td width="98" style="width:73.25pt;border:solid windowtext 1.0pt;mso-border-alt:
  solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:16.15pt">
  <p class="MsoNormal" align="center" style="text-align:center;punctuation-wrap:
  simple;text-autospace:none;mso-line-break-override:restrictions;mso-yfti-cnfc:
  517;mso-element:frame;mso-element-frame-hspace:9.0pt;mso-element-wrap:around;
  mso-element-anchor-vertical:paragraph;mso-element-anchor-horizontal:margin;
  mso-element-top:.05pt;mso-height-rule:exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;color:#1F2328;background:white">Dataset<span style="text-transform:uppercase"><o:p></o:p></span></span></p>
  </td>
  <td width="69" style="width:51.4pt;border:solid windowtext 1.0pt;border-left:
  none;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:16.15pt">
  <p class="MsoNormal" align="center" style="text-align:center;punctuation-wrap:
  simple;text-autospace:none;mso-line-break-override:restrictions;mso-yfti-cnfc:
  1;mso-element:frame;mso-element-frame-hspace:9.0pt;mso-element-wrap:around;
  mso-element-anchor-vertical:paragraph;mso-element-anchor-horizontal:margin;
  mso-element-top:.05pt;mso-height-rule:exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;color:#1F2328;background:white">Instance<span style="text-transform:uppercase"><o:p></o:p></span></span></p>
  </td>
  <td width="66" style="width:49.55pt;border:solid windowtext 1.0pt;border-left:
  none;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:16.15pt">
  <p class="MsoNormal" align="center" style="text-align:center;punctuation-wrap:
  simple;text-autospace:none;mso-line-break-override:restrictions;mso-yfti-cnfc:
  1;mso-element:frame;mso-element-frame-hspace:9.0pt;mso-element-wrap:around;
  mso-element-anchor-vertical:paragraph;mso-element-anchor-horizontal:margin;
  mso-element-top:.05pt;mso-height-rule:exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;color:#1F2328;background:white">Clusters<span style="text-transform:uppercase"><o:p></o:p></span></span></p>
  </td>
  <td width="128" style="width:96.0pt;border:solid windowtext 1.0pt;border-left:
  none;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:16.15pt">
  <p class="MsoNormal" align="center" style="text-align:center;punctuation-wrap:
  simple;text-autospace:none;mso-line-break-override:restrictions;mso-yfti-cnfc:
  1;mso-element:frame;mso-element-frame-hspace:9.0pt;mso-element-wrap:around;
  mso-element-anchor-vertical:paragraph;mso-element-anchor-horizontal:margin;
  mso-element-top:.05pt;mso-height-rule:exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;color:#1F2328;background:white">Views(dimension)<span style="text-transform:uppercase"><o:p></o:p></span></span></p>
  </td>
  <td width="174" style="width:130.25pt;border:solid windowtext 1.0pt;border-left:
  none;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:16.15pt">
  <p class="MsoNormal" align="center" style="text-align:center;punctuation-wrap:
  simple;text-autospace:none;mso-line-break-override:restrictions;mso-yfti-cnfc:
  1;mso-element:frame;mso-element-frame-hspace:9.0pt;mso-element-wrap:around;
  mso-element-anchor-vertical:paragraph;mso-element-anchor-horizontal:margin;
  mso-element-top:.05pt;mso-height-rule:exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;color:#1F2328;background:white">Description<span style="text-transform:uppercase"><o:p></o:p></span></span></p>
  </td>
  <td width="55" style="width:41.55pt;border:solid windowtext 1.0pt;border-left:
  none;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:16.15pt">
  <p class="MsoNormal" align="center" style="text-align:center;punctuation-wrap:
  simple;text-autospace:none;mso-line-break-override:restrictions;mso-yfti-cnfc:
  1;mso-element:frame;mso-element-frame-hspace:9.0pt;mso-element-wrap:around;
  mso-element-anchor-vertical:paragraph;mso-element-anchor-horizontal:margin;
  mso-element-top:.05pt;mso-height-rule:exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;color:#1F2328;background:white">Type<span style="text-transform:uppercase"><o:p></o:p></span></span></p>
  </td>
 </tr>
 <tr style="mso-yfti-irow:0;height:118.45pt">
  <td width="98" style="width:73.25pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  background:#F2F2F2;mso-background-themecolor:background1;mso-background-themeshade:
  242;padding:0cm 5.4pt 0cm 5.4pt;height:118.45pt">
  <p class="MsoNormal" align="center" style="text-align:center;punctuation-wrap:
  simple;text-autospace:none;mso-line-break-override:restrictions;mso-yfti-cnfc:
  68;mso-element:frame;mso-element-frame-hspace:9.0pt;mso-element-wrap:around;
  mso-element-anchor-vertical:paragraph;mso-element-anchor-horizontal:margin;
  mso-element-top:.05pt;mso-height-rule:exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;color:#1F2328">Yale<b><span style="text-transform:uppercase;background:white"><o:p></o:p></span></b></span></p>
  </td>
  <td width="69" style="width:51.4pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;background:#F2F2F2;mso-background-themecolor:background1;
  mso-background-themeshade:242;padding:0cm 5.4pt 0cm 5.4pt;height:118.45pt">
  <p class="MsoNormal" style="punctuation-wrap:simple;text-autospace:none;
  mso-line-break-override:restrictions;mso-yfti-cnfc:64;mso-element:frame;
  mso-element-frame-hspace:9.0pt;mso-element-wrap:around;mso-element-anchor-vertical:
  paragraph;mso-element-anchor-horizontal:margin;mso-element-top:.05pt;
  mso-height-rule:exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;
  color:#1F2328">165<o:p></o:p></span></p>
  </td>
  <td width="66" style="width:49.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;background:#F2F2F2;mso-background-themecolor:
  background1;mso-background-themeshade:242;padding:0cm 5.4pt 0cm 5.4pt;
  height:118.45pt">
  <p class="MsoNormal" style="punctuation-wrap:simple;text-autospace:none;
  mso-line-break-override:restrictions;mso-yfti-cnfc:64;mso-element:frame;
  mso-element-frame-hspace:9.0pt;mso-element-wrap:around;mso-element-anchor-vertical:
  paragraph;mso-element-anchor-horizontal:margin;mso-element-top:.05pt;
  mso-height-rule:exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;
  color:#1F2328">15<o:p></o:p></span></p>
  </td>
  <td width="128" style="width:96.0pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;background:#F2F2F2;mso-background-themecolor:
  background1;mso-background-themeshade:242;padding:0cm 5.4pt 0cm 5.4pt;
  height:118.45pt">
  <p class="MsoNormal" style="punctuation-wrap:simple;text-autospace:none;
  mso-line-break-override:restrictions;mso-yfti-cnfc:64;mso-element:frame;
  mso-element-frame-hspace:9.0pt;mso-element-wrap:around;mso-element-anchor-vertical:
  paragraph;mso-element-anchor-horizontal:margin;mso-element-top:.05pt;
  mso-height-rule:exactly"><span class="GramE"><span lang="EN-US" style="font-family:
  &quot;Calibri&quot;,sans-serif;color:#1F2328">Intensity(</span></span><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;color:#1F2328">4096), <span class="GramE">LBP(</span>3304), <span class="GramE">Gabor(</span>6750)<o:p></o:p></span></p>
  </td>
  <td width="174" style="width:130.25pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;background:#F2F2F2;mso-background-themecolor:
  background1;mso-background-themeshade:242;padding:0cm 5.4pt 0cm 5.4pt;
  height:118.45pt">
  <p class="MsoNormal" align="left" style="text-align:left;page-break-before:always;
  mso-pagination:widow-orphan lines-together;page-break-after:avoid;punctuation-trim:
  leading;punctuation-wrap:simple;vertical-align:middle;mso-yfti-cnfc:64;
  mso-element:frame;mso-element-frame-hspace:9.0pt;mso-element-wrap:around;
  mso-element-anchor-vertical:paragraph;mso-element-anchor-horizontal:margin;
  mso-element-top:.05pt;mso-height-rule:exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;color:#1F2328">This dataset consists
  of 165&nbsp;gray-scale face images belonging to 15 subjects with each subject
  containing 11 images [1].<o:p></o:p></span></p>
  </td>
  <td width="55" style="width:41.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;background:#F2F2F2;mso-background-themecolor:
  background1;mso-background-themeshade:242;padding:0cm 5.4pt 0cm 5.4pt;
  height:118.45pt">
  <p class="MsoNormal" style="punctuation-wrap:simple;text-autospace:none;
  mso-line-break-override:restrictions;mso-yfti-cnfc:64;mso-element:frame;
  mso-element-frame-hspace:9.0pt;mso-element-wrap:around;mso-element-anchor-vertical:
  paragraph;mso-element-anchor-horizontal:margin;mso-element-top:.05pt;
  mso-height-rule:exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;
  color:#1F2328">Face<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style="mso-yfti-irow:1;height:115.6pt">
  <td width="98" style="width:73.25pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:115.6pt">
  <p class="MsoNormal" align="center" style="text-align:center;punctuation-wrap:
  simple;text-autospace:none;mso-line-break-override:restrictions;mso-yfti-cnfc:
  4;mso-element:frame;mso-element-frame-hspace:9.0pt;mso-element-wrap:around;
  mso-element-anchor-vertical:paragraph;mso-element-anchor-horizontal:margin;
  mso-element-top:.05pt;mso-height-rule:exactly"><span class="SpellE"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;color:#1F2328;background:
  white">HandWritten</span></span><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;
  color:#1F2328;background:white"><o:p></o:p></span></p>
  </td>
  <td width="69" style="width:51.4pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:115.6pt">
  <p class="MsoNormal" style="punctuation-wrap:simple;text-autospace:none;
  mso-line-break-override:restrictions;mso-element:frame;mso-element-frame-hspace:
  9.0pt;mso-element-wrap:around;mso-element-anchor-vertical:paragraph;
  mso-element-anchor-horizontal:margin;mso-element-top:.05pt;mso-height-rule:
  exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;color:#1F2328;
  background:white">2000<o:p></o:p></span></p>
  </td>
  <td width="66" style="width:49.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:115.6pt">
  <p class="MsoNormal" style="punctuation-wrap:simple;text-autospace:none;
  mso-line-break-override:restrictions;mso-element:frame;mso-element-frame-hspace:
  9.0pt;mso-element-wrap:around;mso-element-anchor-vertical:paragraph;
  mso-element-anchor-horizontal:margin;mso-element-top:.05pt;mso-height-rule:
  exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;color:#1F2328;
  background:white">10<o:p></o:p></span></p>
  </td>
  <td width="128" style="width:96.0pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:115.6pt">
  <p class="MsoNormal" style="punctuation-wrap:simple;text-autospace:none;
  mso-line-break-override:restrictions;mso-element:frame;mso-element-frame-hspace:
  9.0pt;mso-element-wrap:around;mso-element-anchor-vertical:paragraph;
  mso-element-anchor-horizontal:margin;mso-element-top:.05pt;mso-height-rule:
  exactly"><span class="GramE"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;
  color:#1F2328;background:white">FOU(</span></span><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;color:#1F2328;background:white">76), <span class="GramE">FAC(</span>216), <span class="GramE">KAR(</span>64), <span class="GramE">PIX(</span>240), <span class="GramE">ZER(</span>47), <span class="GramE">MOR(</span>6)<o:p></o:p></span></p>
  </td>
  <td width="174" style="width:130.25pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:115.6pt">
  <p class="MsoNormal" style="punctuation-wrap:simple;text-autospace:none;
  mso-line-break-override:restrictions;mso-element:frame;mso-element-frame-hspace:
  9.0pt;mso-element-wrap:around;mso-element-anchor-vertical:paragraph;
  mso-element-anchor-horizontal:margin;mso-element-top:.05pt;mso-height-rule:
  exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;color:#1F2328;
  background:white">This dataset consists of features of handwritten numerals
  ('0'--'9') extracted from a collection of Dutch utility maps [2].<o:p></o:p></span></p>
  </td>
  <td width="55" style="width:41.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:115.6pt">
  <p class="MsoNormal" style="punctuation-wrap:simple;text-autospace:none;
  mso-line-break-override:restrictions;mso-element:frame;mso-element-frame-hspace:
  9.0pt;mso-element-wrap:around;mso-element-anchor-vertical:paragraph;
  mso-element-anchor-horizontal:margin;mso-element-top:.05pt;mso-height-rule:
  exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;color:#1F2328;
  background:white">Image<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style="mso-yfti-irow:2;height:108.55pt">
  <td width="98" style="width:73.25pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  background:#F2F2F2;mso-background-themecolor:background1;mso-background-themeshade:
  242;padding:0cm 5.4pt 0cm 5.4pt;height:108.55pt">
  <p class="MsoNormal" align="center" style="text-align:center;punctuation-wrap:
  simple;text-autospace:none;mso-line-break-override:restrictions;mso-yfti-cnfc:
  68;mso-element:frame;mso-element-frame-hspace:9.0pt;mso-element-wrap:around;
  mso-element-anchor-vertical:paragraph;mso-element-anchor-horizontal:margin;
  mso-element-top:.05pt;mso-height-rule:exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;color:black;mso-color-alt:windowtext">3Sources</span><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;color:#1F2328"><o:p></o:p></span></p>
  </td>
  <td width="69" style="width:51.4pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;background:#F2F2F2;mso-background-themecolor:background1;
  mso-background-themeshade:242;padding:0cm 5.4pt 0cm 5.4pt;height:108.55pt">
  <p class="MsoNormal" style="punctuation-wrap:simple;text-autospace:none;
  mso-line-break-override:restrictions;mso-yfti-cnfc:64;mso-element:frame;
  mso-element-frame-hspace:9.0pt;mso-element-wrap:around;mso-element-anchor-vertical:
  paragraph;mso-element-anchor-horizontal:margin;mso-element-top:.05pt;
  mso-height-rule:exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;
  color:#1F2328">169<o:p></o:p></span></p>
  </td>
  <td width="66" style="width:49.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;background:#F2F2F2;mso-background-themecolor:
  background1;mso-background-themeshade:242;padding:0cm 5.4pt 0cm 5.4pt;
  height:108.55pt">
  <p class="MsoNormal" style="punctuation-wrap:simple;text-autospace:none;
  mso-line-break-override:restrictions;mso-yfti-cnfc:64;mso-element:frame;
  mso-element-frame-hspace:9.0pt;mso-element-wrap:around;mso-element-anchor-vertical:
  paragraph;mso-element-anchor-horizontal:margin;mso-element-top:.05pt;
  mso-height-rule:exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;
  color:#1F2328">6<o:p></o:p></span></p>
  </td>
  <td width="128" style="width:96.0pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;background:#F2F2F2;mso-background-themecolor:
  background1;mso-background-themeshade:242;padding:0cm 5.4pt 0cm 5.4pt;
  height:108.55pt">
  <p class="MsoNormal" style="punctuation-wrap:simple;text-autospace:none;
  mso-line-break-override:restrictions;mso-yfti-cnfc:64;mso-element:frame;
  mso-element-frame-hspace:9.0pt;mso-element-wrap:around;mso-element-anchor-vertical:
  paragraph;mso-element-anchor-horizontal:margin;mso-element-top:.05pt;
  mso-height-rule:exactly"><span class="GramE"><span lang="EN-US" style="font-family:
  &quot;Calibri&quot;,sans-serif;color:#1F2328">Reuters(</span></span><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;color:#1F2328">3068), <span class="GramE">BBC(</span>3560), <span class="GramE">Guardian(</span>3631)<o:p></o:p></span></p>
  </td>
  <td width="174" style="width:130.25pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;background:#F2F2F2;mso-background-themecolor:
  background1;mso-background-themeshade:242;padding:0cm 5.4pt 0cm 5.4pt;
  height:108.55pt">
  <p class="MsoNormal" style="punctuation-wrap:simple;text-autospace:none;
  mso-line-break-override:restrictions;mso-yfti-cnfc:64;mso-element:frame;
  mso-element-frame-hspace:9.0pt;mso-element-wrap:around;mso-element-anchor-vertical:
  paragraph;mso-element-anchor-horizontal:margin;mso-element-top:.05pt;
  mso-height-rule:exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;
  color:#1F2328">A new multi-view text dataset collected from three well-known
  online news sources: BBC, Reuters, and The Guardian [3].<o:p></o:p></span></p>
  </td>
  <td width="55" style="width:41.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;background:#F2F2F2;mso-background-themecolor:
  background1;mso-background-themeshade:242;padding:0cm 5.4pt 0cm 5.4pt;
  height:108.55pt">
  <p class="MsoNormal" style="punctuation-wrap:simple;text-autospace:none;
  mso-line-break-override:restrictions;mso-yfti-cnfc:64;mso-element:frame;
  mso-element-frame-hspace:9.0pt;mso-element-wrap:around;mso-element-anchor-vertical:
  paragraph;mso-element-anchor-horizontal:margin;mso-element-top:.05pt;
  mso-height-rule:exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;
  color:#1F2328">Text<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style="mso-yfti-irow:3;height:143.65pt">
  <td width="98" style="width:73.25pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0cm 5.4pt 0cm 5.4pt;height:143.65pt">
  <p class="MsoNormal" align="center" style="text-align:center;punctuation-wrap:
  simple;text-autospace:none;mso-line-break-override:restrictions;mso-yfti-cnfc:
  4;mso-element:frame;mso-element-frame-hspace:9.0pt;mso-element-wrap:around;
  mso-element-anchor-vertical:paragraph;mso-element-anchor-horizontal:margin;
  mso-element-top:.05pt;mso-height-rule:exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif">Cora<span style="color:#1F2328"><o:p></o:p></span></span></p>
  </td>
  <td width="69" style="width:51.4pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:143.65pt">
  <p class="MsoNormal" style="punctuation-wrap:simple;text-autospace:none;
  mso-line-break-override:restrictions;mso-element:frame;mso-element-frame-hspace:
  9.0pt;mso-element-wrap:around;mso-element-anchor-vertical:paragraph;
  mso-element-anchor-horizontal:margin;mso-element-top:.05pt;mso-height-rule:
  exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;color:#1F2328">2708<o:p></o:p></span></p>
  </td>
  <td width="66" style="width:49.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:143.65pt">
  <p class="MsoNormal" style="punctuation-wrap:simple;text-autospace:none;
  mso-line-break-override:restrictions;mso-element:frame;mso-element-frame-hspace:
  9.0pt;mso-element-wrap:around;mso-element-anchor-vertical:paragraph;
  mso-element-anchor-horizontal:margin;mso-element-top:.05pt;mso-height-rule:
  exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;color:#1F2328">7<o:p></o:p></span></p>
  </td>
  <td width="128" style="width:96.0pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:143.65pt">
  <p class="MsoNormal" style="punctuation-wrap:simple;text-autospace:none;
  mso-line-break-override:restrictions;mso-element:frame;mso-element-frame-hspace:
  9.0pt;mso-element-wrap:around;mso-element-anchor-vertical:paragraph;
  mso-element-anchor-horizontal:margin;mso-element-top:.05pt;mso-height-rule:
  exactly"><span class="GramE"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;
  color:#1F2328">Content(</span></span><span lang="EN-US" style="font-family:
  &quot;Calibri&quot;,sans-serif;color:#1F2328">1433), <span class="GramE">Inbound(</span>2708),
  <span class="GramE">Outbound(</span>2708), <span class="GramE">Cites(</span>2708)<o:p></o:p></span></p>
  </td>
  <td width="174" style="width:130.25pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:143.65pt">
  <p class="MsoNormal" style="punctuation-wrap:simple;text-autospace:none;
  mso-line-break-override:restrictions;mso-element:frame;mso-element-frame-hspace:
  9.0pt;mso-element-wrap:around;mso-element-anchor-vertical:paragraph;
  mso-element-anchor-horizontal:margin;mso-element-top:.05pt;mso-height-rule:
  exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;color:#1F2328">The
  archive contains 2708 documents over the 7 labels (<span class="SpellE">Neural_Networks</span>,
  <span class="SpellE">Rule_Learning</span>, <span class="SpellE">Reinforcement_Learning</span>,
  <span class="SpellE">Probabilistic_Methods</span>, Theory, <span class="SpellE">Genetic_Algorithms</span>,
  <span class="SpellE">Case_Based</span>) [4].<o:p></o:p></span></p>
  </td>
  <td width="55" style="width:41.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0cm 5.4pt 0cm 5.4pt;height:143.65pt">
  <p class="MsoNormal" style="punctuation-wrap:simple;text-autospace:none;
  mso-line-break-override:restrictions;mso-element:frame;mso-element-frame-hspace:
  9.0pt;mso-element-wrap:around;mso-element-anchor-vertical:paragraph;
  mso-element-anchor-horizontal:margin;mso-element-top:.05pt;mso-height-rule:
  exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;color:#1F2328">Text<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style="mso-yfti-irow:4;mso-yfti-lastrow:yes;height:73.95pt">
  <td width="98" style="width:73.25pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  background:#F2F2F2;mso-background-themecolor:background1;mso-background-themeshade:
  242;padding:0cm 5.4pt 0cm 5.4pt;height:73.95pt">
  <p class="MsoNormal" align="center" style="text-align:center;punctuation-wrap:
  simple;text-autospace:none;mso-line-break-override:restrictions;mso-yfti-cnfc:
  68;mso-element:frame;mso-element-frame-hspace:9.0pt;mso-element-wrap:around;
  mso-element-anchor-vertical:paragraph;mso-element-anchor-horizontal:margin;
  mso-element-top:.05pt;mso-height-rule:exactly"><span class="SpellE"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;color:black;mso-color-alt:
  windowtext">Citeseer</span></span><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;
  color:#1F2328;text-transform:uppercase"><o:p></o:p></span></p>
  </td>
  <td width="69" style="width:51.4pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;background:#F2F2F2;mso-background-themecolor:background1;
  mso-background-themeshade:242;padding:0cm 5.4pt 0cm 5.4pt;height:73.95pt">
  <p class="MsoNormal" style="punctuation-wrap:simple;text-autospace:none;
  mso-line-break-override:restrictions;mso-yfti-cnfc:64;mso-element:frame;
  mso-element-frame-hspace:9.0pt;mso-element-wrap:around;mso-element-anchor-vertical:
  paragraph;mso-element-anchor-horizontal:margin;mso-element-top:.05pt;
  mso-height-rule:exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;
  color:#1F2328">3312<o:p></o:p></span></p>
  </td>
  <td width="66" style="width:49.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;background:#F2F2F2;mso-background-themecolor:
  background1;mso-background-themeshade:242;padding:0cm 5.4pt 0cm 5.4pt;
  height:73.95pt">
  <p class="MsoNormal" style="punctuation-wrap:simple;text-autospace:none;
  mso-line-break-override:restrictions;mso-yfti-cnfc:64;mso-element:frame;
  mso-element-frame-hspace:9.0pt;mso-element-wrap:around;mso-element-anchor-vertical:
  paragraph;mso-element-anchor-horizontal:margin;mso-element-top:.05pt;
  mso-height-rule:exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;
  color:#1F2328">6<o:p></o:p></span></p>
  </td>
  <td width="128" style="width:96.0pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;background:#F2F2F2;mso-background-themecolor:
  background1;mso-background-themeshade:242;padding:0cm 5.4pt 0cm 5.4pt;
  height:73.95pt">
  <p class="MsoNormal" style="punctuation-wrap:simple;text-autospace:none;
  mso-line-break-override:restrictions;mso-yfti-cnfc:64;mso-element:frame;
  mso-element-frame-hspace:9.0pt;mso-element-wrap:around;mso-element-anchor-vertical:
  paragraph;mso-element-anchor-horizontal:margin;mso-element-top:.05pt;
  mso-height-rule:exactly"><span class="GramE"><span lang="EN-US" style="font-family:
  &quot;Calibri&quot;,sans-serif;color:#1F2328">Content(</span></span><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;color:#1F2328">3703), <span class="GramE">Cites(</span>4732)<o:p></o:p></span></p>
  </td>
  <td width="174" style="width:130.25pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;background:#F2F2F2;mso-background-themecolor:
  background1;mso-background-themeshade:242;padding:0cm 5.4pt 0cm 5.4pt;
  height:73.95pt">
  <p class="MsoNormal" style="punctuation-wrap:simple;text-autospace:none;
  mso-line-break-override:restrictions;mso-yfti-cnfc:64;mso-element:frame;
  mso-element-frame-hspace:9.0pt;mso-element-wrap:around;mso-element-anchor-vertical:
  paragraph;mso-element-anchor-horizontal:margin;mso-element-top:.05pt;
  mso-height-rule:exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;
  color:#1F2328">The archive contains 3312 documents over the 6 labels (Agents,
  IR, DB, AI, HCI, ML) [4]<o:p></o:p></span></p>
  </td>
  <td width="55" style="width:41.55pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;background:#F2F2F2;mso-background-themecolor:
  background1;mso-background-themeshade:242;padding:0cm 5.4pt 0cm 5.4pt;
  height:73.95pt">
  <p class="MsoNormal" style="punctuation-wrap:simple;text-autospace:none;
  mso-line-break-override:restrictions;mso-yfti-cnfc:64;mso-element:frame;
  mso-element-frame-hspace:9.0pt;mso-element-wrap:around;mso-element-anchor-vertical:
  paragraph;mso-element-anchor-horizontal:margin;mso-element-top:.05pt;
  mso-height-rule:exactly"><span lang="EN-US" style="font-family:&quot;Calibri&quot;,sans-serif;
  color:#1F2328">Text<o:p></o:p></span></p>
  </td>
 </tr>
</tbody></table>

[1] M.-S. Chen, C.-D. Wang, and J.-H. Lai, “Low-rank Tensor Based Proximity Learning for Multi-view Clustering,” IEEE Transactions on Knowledge and Data Engineering, pp. 1–1, Jan. 2022, doi: 10.1109/TKDE.2022.3151861.
[2] F. Nie, L. Tian, and X. Li, “Multiview clustering via adaptively weighted procrustes,” in Proc. ACM Int. Conf. Knowl. Discov. Data Min., 2018, pp. 2022–2030. doi: 10.1145/3219819.3220049.
[3] H. Wei, L. Chen, C. L. P. Chen, J. Duan, R. Han, and L. Guo, “Fuzzy clustering for multiview data by combining latent information,” Applied Soft Computing, vol. 126, p. 109140, Sep. 2022, doi: 10.1016/j.asoc.2022.109140.
[4] S.-G. Fang, D. Huang, X.-S. Cai, C.-D. Wang, C. He, and Y. Tang, “Efficient Multi-view Clustering via Unified and Discrete Bipartite Graph Learning,” IEEE Transactions on Neural Networks and Learning Systems, pp. 1–12, Apr. 2023, doi: 10.1109/TNNLS.2023.3261460.
