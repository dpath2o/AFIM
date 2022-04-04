
# Table of Contents

1.  [Background](#orgba38d4f):personal:
2.  [Administration](#orgce534f4)
3.  [Research](#orgded9e8d):res:
    1.  [UTas Library Services](#org097d2c2):bookmarks:
    2.  [Modelling](#org2bb8890):modelling:
        1.  [CICE](#orgfb12bc2):coupled:seaice:modelling:
        2.  [WAVEWATCH III](#orge8bceea):bookmark:modelling:ocean:
        3.  [SWAN (Simulating WAves Nearshore)](#orgc7ca88c):bookmark:modelling:ocean:
        4.  [Weather Research & Forecasting (WRF) Model](#orgd8641cf):bookmark:modelling:atmo:
        5.  [Model Coupling Toolkit (MCT)](#orga32fc3e):bookmark:modelling:atmo:
        6.  [COAWST](#orgfaec473):modelling:ocean:coupled:
        7.  [ROMS](#org926b2c9):ocean:modelling:
    3.  [Data](#org8dc835f):data:
    4.  [Research Plan 2.0](#org21e5ac1):plan:
    5.  [Candidature](#org781c4c5):report:
        1.  [Meetings](#orgf8326cd)
        2.  [Ethics Approval](#org351a2f5)
        3.  [Confirmation of Candidature](#org763433d)
        4.  [Annual Review of Candidature](#org07ffddb)
        5.  [Thesis Submission and Examination](#org76c432d)
        6.  [Step 1: Abstract](#org3907061)
        7.  [Step 2: Thesis/Exegesis Submission](#org443df34)
        8.  [Step 3: Examiner Reports](#org0745770)
        9.  [Step 4: Final Thesis/Exegesis Submission](#orgd060575)
        10. [Step 5: Graduation](#orgbf8fed8)
4.  [Meetings and Correspondence](#org3a8f892):mtg:email:
5.  [Computing](#orgf5ae7f9):comp:

<span class="underline">Technical Note</span>: This file is written in [Emacs-Org](https://orgmode.org) (which, as an
aside, I have been enamoured with for more years than I care to
admit). As such it is exported to the [markdown](https://www.markdownguide.org) language that is
favoured by <https://github.com>. However if one would like a slightly
more native experience to Emacs-Org then go [here](./README.md). Irrespective, this
README.md will be clobbered with each edit of its parent org file and
so if you are wanting to edit, edit the org file :) 


<a id="orgba38d4f"></a>

# Background     :personal:

At the start of 2020 I decided to pursue a PhD in oceanography that
comprised [ocean modelling](http://www.cmar.csiro.au/staff/oke/pubs/England_and_Oke_2001.pdf), the [Southern Ocean](https://tos.org/oceanography/issue/volume-25-issue-03) and [Antarctica](https://www.scar.org). This was
spurned largely from my [previous](https://www.cencoos.org) [professional life](http://imos.org.au)
[coastal oceanographer](https://scripps.ucsd.edu/research/topics/coastal-oceanography). That previous pursuit centred mainly around a remote
sensing technology called [high frequency radar](https://tos.org/oceanography/assets/docs/10-2_paduan1.pdf) which is, and I used, in a
broad number of ways to understand the dynamics (the motion of the
of the upper ocean) through the digital application of signal
processing of the Physical representaiton of a Doppler shifted Bragg
frequency. After roughly 15 years in that field I left to diversify my
career skillset while ticking some service-related philosophies about
being [contributing to a new nation](https://en.wikipedia.org/wiki/National_service). After four years of [learning how-to drive Navy ships](https://www.navy.gov.au/sites/default/files/documents/Warfare_Officers_Career_Handbook.pdf) I
specialised as a [Meteorologic and Oceanographic Officer](https://www.defencejobs.gov.au/jobs/reserves/navy/meteorologist-and-oceanographer) in that
[organisation](https://www.navy.gov.au) and [operational-ised](https://www.youtube.com/watch?v=_1roFUwV7ss) my [scientific skillset](https://oceansci.ucsc.edu/academics/graduate/ms.html). At about the
same time I became fascinated with Antarctica and its criticality both
[climatically](https://tos.org/oceanography/article/southern-ocean-warming) and [strategically](https://www.antarctica.gov.au/about-us/antarctic-strategy-and-action-plan/). Hence I began searching for a PhD
project that would fulfil my interest.

Fortunately, it did not take too long before I was introduced to [Alex Fraser](https://tasmanian.com.au/stories/alex-fraser/) and [a project](./ResearchPlan/project_proposal/PROJECT_PROPOSAL.pdf) that he had on his digital shelf that from the outset ticked all my
interest boxes. I can't recall exactly where or when in I was introduced to the various forms of
[sea ice](https://en.wikipedia.org/wiki/Sea_ice) and [its role in polar oceanography](https://tos.org/oceanography/issue/volume-24-issue-03) but I do remember it
being very [Arctic-centric](http://nsidc.org/arcticseaicenews/) with this vague concept that [landfast sea ice (fast ice)](https://arctic.noaa.gov/Report-Card/Report-Card-2018/ArtMID/7878/ArticleID/788/Landfast-Sea-Ice-in-a-Changing-Arctic) played an important role in Arctic sea ice. So when Alex
provided me with a rough outline of the project he wanted to pursue
modelling fast ice and it being [largely neglected thus far](https://www.sciencedirect.com/user/identity/landing?code=2QD--DJX9fNQGcwrhJdFSQcbTnjTstfb_oMPFrXr&state=retryCounter%3D0%26csrfToken%3D282ae9df-d109-4361-a158-7399407fb2cd%26idpPolicy%3Durn%253Acom%253Aelsevier%253Aidp%253Apolicy%253Aproduct%253Ainst_assoc%26returnUrl%3D%252Fscience%252Farticle%252Fpii%252FS1463500321001712%26prompt%3Dnone%26cid%3Darp-05cb79cb-d88b-43bd-92e0-7af2740528af) in
circumpolar modelling efforts, I knew the meat on this bone was likely
well marbled.

So in mid-2021 I was accepted and enrolled as a part-time PhD student at
the [University of Tasmania](https://www.googleadservices.com/pagead/aclk?sa=L&ai=DChcSEwjb1_Cu-_j2AhWFmmYCHS-oCUAYABAAGgJzbQ&ohost=www.google.com&cid=CAESbeD2zyC-aH5iZzoJORucLfQbhrhY_ooy9py4o71TC_kz3xRqKOproNUbUwtZtC5573rDHXbolu_5OSzJ1IBTOd5vLU41bby25TJOXU74Xos06Noqir-L1xomT419-70EUgr0C5r7Znshu44t3TM&sig=AOD64_2ugbp2WvnZ4lnJ32JarlVkJGBBcA&q&adurl&ved=2ahUKEwjJxumu-_j2AhXm63MBHchuCCMQ0Qx6BAgDEAE) [Institute of Marine Science](https://www.imas.utas.edu.au) through the
[Australian Antarctic Partnership Program](https://aappartnership.org.au). I spent the second-half of
that year [reading](./references) and constructing an initial [Reasearch Plan](./ResearchPlan/doc/researchplan.pdf). However, as all good
projects evolve, it became apparent that I should be aligning my
project with a well-supported Australian sea ice modelling effort, not
only for my own support, but also to aim for a more impactful goal of
incorporating fast ice into a nationally recognised/supported
model.

After a brief interlude in the first quarter of 2022 I'm now pursuing
fast ice modelling with an eye towards incorporating it into [COSIMA](http://cosima.org.au).


<a id="orgce534f4"></a>

# Administration

<./admin/admin.md>


<a id="orgded9e8d"></a>

# Research     :res:


<a id="org097d2c2"></a>

## UTas Library Services     :bookmarks:

[UTas Research](https://www.utas.edu.au/library/research)


<a id="org2bb8890"></a>

## Modelling     :modelling:


<a id="orgfb12bc2"></a>

### CICE     :coupled:seaice:modelling:

<./MODELS/cice.md>


<a id="orge8bceea"></a>

### WAVEWATCH III     :bookmark:modelling:ocean:

<./MODELS/wwIII.md>


<a id="orgc7ca88c"></a>

### SWAN (Simulating WAves Nearshore)     :bookmark:modelling:ocean:

<./MODELS/swan.md>


<a id="orgd8641cf"></a>

### Weather Research & Forecasting (WRF) Model     :bookmark:modelling:atmo:

<./MODELS/wrf.md>


<a id="orga32fc3e"></a>

### Model Coupling Toolkit (MCT)     :bookmark:modelling:atmo:

<./MODELS/mct.md>


<a id="orgfaec473"></a>

### COAWST     :modelling:ocean:coupled:

<./MODELS/coawst.md>


<a id="org926b2c9"></a>

### ROMS     :ocean:modelling:

<./MODELS/metroms.md>


<a id="org8dc835f"></a>

## Data     :data:

<./data.md>


<a id="org21e5ac1"></a>

## Research Plan 2.0     :plan:

Being a dynamic document this simply a new version of the previous
plan. As it should it will be kept separate from this document and can
be found [here](./researchplan/researchplan.md) (as an Emacs-Org)


<a id="org781c4c5"></a>

## Candidature     :report:


<a id="orgf8326cd"></a>

### Meetings

You are encouraged to meet with your supervisory team on a regular basis (e.g. fortnightly) to discuss project progress, troubleshoot issues and develop action plans for future research. These regular meetings support fruitful communication between you and your team, and assists all parties to manage their expectations.
Meeting Agendas: You are responsible for preparing the agenda for your meetings and documenting the discussion and outcomes. This record provides clarity for all parties. The agenda can include:
Meeting date/attendees

-   Progress on matters arising from the last meeting
-   Other project progress since the last meeting
-   Matters for discussion/troubleshooting
-   Progress towards your next candidature milestone(s)
-   Actions/outcomes arising from current meeting
-   Any other business


<a id="org351a2f5"></a>

### Ethics Approval

The University of Tasmania (UTAS) is dedicated to creating and maintaining an environment that promotes the responsible and ethical conduct of research in accordance with the Australian Code for the Responsible Conduct of Research.
It is essential that you discuss with your supervisors at the very outset whether ethics approval is needed and how to obtain it. Many research projects, including Doctor of Philosophy and Master of Research projects, require Human Ethics or Animal Ethics approval to proceed.
The Research Integrity and Ethics Unit is there to support researchers at UTAS. Please contact them with any questions relating to an ethics application, a clinical trial or a research integrity issue.


<a id="org763433d"></a>

### Confirmation of Candidature

The purpose of Confirmation of Candidature is to assess:

-   Your academic preparedness and whether you have developed a clearly defined,
    coherent and feasible research project and whether you have documented this
    adequately in your research plan
-   Whether you have met coursework requirements as specified in your letter of offer
-   Whether other specific requirements (e.g. obtaining ethics approval) have been met
-   The suitability of the supervisory team to support you to completion
-   Whether your oral and written skills are sufficient for undertaking a Higher Degree by Research
-   The likelihood of you completing your Higher Degree by Research within the maximum degree period
-   Discuss the requirements of your Confirmation of Candidature with your
    supervisors and Graduate Research Coordinator (GRC) as early as possible

In order to have your candidature confirmed, you must fulfil the following
requirements:

-   Research Plan
-   Ethics Requirements
-   Data Management Plan
-   Written Work
-   Oral Presentation
-   Peer Review
-   Supervisor-Candidate meetings
-   Coursework
-   Any additional criteria required by your school or specified in your letter of offer

Detailed information can be found in the Research Training Ordinance and the
Interim HDR Candidature Management procedure and Interim HDR Reviews of Progress
procedure.


<a id="org07ffddb"></a>

### Annual Review of Candidature

Following successful Confirmation of Candidature, your progress will be formally
reviewed every 12 calendar months through candidature until submission.
The purpose of the annual review is to make an assessment of:

-   Your academic performance consistent with your research plan
-   The adequacy of research infrastructure and resources (including the
    supervisory team relationship) needed to complete your research project within
    the maximum degree period

You will be required to complete the Annual Review tab in iGRad. This allows you
to fill out information and upload documents regarding your Annual Review,
before a meeting is held with your Graduate Research Coordinator (GRC).
Detailed information can be found in the Research Training Ordinance and the
Interim HDR Candidature Management procedure and Interim HDR Reviews of Progress
procedure.


<a id="org76c432d"></a>

### Thesis Submission and Examination


<a id="org3907061"></a>

### Step 1: Abstract

At least eight weeks (8) prior to submitting your thesis/exegesis for
examination, you will need to upload your abstract by completing the Examination
Intention to Submit tab in iGRad.
On receipt of your abstract, a request will be sent to your supervisory team and
Head of School, from the Graduate Research Office, asking them to nominate
examiners for your thesis/exegesis.
The nomination of examiners process is confidential. You are permitted to
provide the names of examiners that you do not want to examine your work (and
accompanying reasons), but that is the extent to which you are involved.
A Chair of Examiners will also be appointed by the Head of School to act as an
independent Chair throughout the examination process. Where the Head of School
is also a member of your supervisory team, the relevant College Dean or
Institute Director will appoint the Chair of Examiners.


<a id="org443df34"></a>

### Step 2: Thesis/Exegesis Submission

Before the maximum expiry date of your candidature, you will need to submit your
thesis/exegesis (and supplementary material if relevant) for examination by
completing the Examination > Thesis Submission tab in iGRad.
If your thesis is less than 60MB in size, and you are only submitting a thesis,
please upload your thesis to the Thesis Submission tab in iGRad.
If your thesis is greater than 60MB in size, or you have supplementary material,
please generate a link to your thesis (and supplementary material) via CloudStor
and then embed the URL link generated via CloudStor onto your title page. Once
completed, please upload the title page (only) to the Thesis Submission tab in
iGRad. To access CloudStor, you can use your UTAS login.
Supplementary material may include, but not be limited to, appendices, research
publications, books, a portfolio of audio/video musical recordings,
documentaries or art works e.g., photos of exhibited items produced during your
work. Please do not duplicate copies of the supplementary material that is
already incorporated into your thesis material.
The thesis/exegesis must include all relevant and signed statements and
declarations before uploading. For further information on these, and other
general requirements, please access the General Requirements of a Higher Degree
by Research Thesis document. Specific details about formatting, length and
referencing will also be dependent upon the discipline in which you have
conducted your research. Please consult with your supervisory team about these
details.
Printed copies of your thesis/exegesis are no longer required, unless
specifically requested by an examiner. In this case, the Graduate Research
Office will request these from you via the School.
Upon receipt of your thesis/exegesis notice will be sent to your Head of School
requesting their approval of your submission.
You will be able to track the examination process via the Examination tab in
iGRad. This will inform you of when examiners have been nominated, when your
thesis/exegesis has been sent and the due dates for examiner reports.


<a id="org0745770"></a>

### Step 3: Examiner Reports

Examiners are requested to return their reports to the Graduate Research Office within six weeks of receiving the thesis/exegesis, however this often takes longer due to workload and personal circumstances. The Graduate Research Office has procedures in place to follow-up on overdue examiner reports.

Once all of the examiner reports have been received, they are sent to the Chair of Examiners. The Chair of Examiners are responsible for providing a recommendation to the Graduate Research Office on the reports, who will then submit their comments and recommendation on the examination outcome to the Dean of Graduate Research.

The Dean of Graduate Research shall consider the recommendations of the examiners and the Chair of Examiners before determining how to proceed to the next stage.

Once determined, we will inform you of the outcome of your examination and send you the examiner reports.

When you have completed any required revisions, please submit your corrected thesis to the Chair of Examiners. They will review your final thesis and, if approved, recommend to the Graduate Research Office that your degree be awarded.

After the Graduate Research Office receive this approval, we can then advise you of the identity of your examiners and invite you to upload your final thesis and graduation requirements by completing the Examination > Graduation Requirements tab in iGRad.


<a id="orgd060575"></a>

### Step 4: Final Thesis/Exegesis Submission

After the Graduate Research Office have invited you to upload your final thesis and graduation requirements, you will need to submit your final thesis/exegesis (including all supplementary material if relevant) by completing the Examination > Graduation Requirements tab in iGRad.

If your thesis is less than 60MB in size, and you are only submitting you’re a thesis, please upload your thesis to the Graduation Requirements tab in iGRad.

If your thesis is greater than 60MB in size, or you have supplementary material, please generate a link to your thesis (and supplementary material) via CloudStor and then embed the URL link generated via CloudStor onto your title page. Once completed, please upload the title page (only) to the Graduation Requirements tab in iGRad. To access CloudStor, you can use your UTAS login.

The thesis/exegesis must include all relevant and signed statements and declarations before uploading. Please be advised that printed copies are no longer required.

To fulfil our graduation requirements, you will need to provide:

Your Final Thesis
All supplementary material (which may include, but not be limited to, appendices, research publications, books, a portfolio of audio/video musical recordings, documentaries or art works e.g., photos of exhibited items produced during your work).
A Thesis Access Form
A thesis abstract for the purpose of AHEGS. This abstract must be strictly 2000 characters (including spaces) or less: This is required for the Australian Higher Education Graduation Statement (AHEGS), a Commonwealth requirement.
For PhD or Professional Doctorate candidates only: A 40-word summary of your thesis to be read at the ceremony to a general audience. Please refer to the 40 Word Summary Guidelines.
Once you have fulfilled the above requirements, the Graduate Research Office will request Academic Senate approval to award your degree. We will advise you when this approval has been received.


<a id="orgbf8fed8"></a>

### Step 5: Graduation


<a id="org3a8f892"></a>

# Meetings and Correspondence     :mtg:email:

<./admin/meetings.md>


<a id="orgf5ae7f9"></a>

# Computing     :comp:

<./admin/computing.md>

