#include <stdio.h>
#include <string.h>
#include "sac.h"

typedef struct {
	int yr;	/* year		*/
	int mo;	/* month	*/
	int day;	/* day		*/
	int hr;	/* hour		*/
	int mn;	/* minute	*/
	float sec;	/* second	*/
} Time;
typedef int bool ; 
bool isleap(int);

static int      days_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31};

#define mod(a,b)	a - ((int)(a/b)) * b

/* fill the time structure *t with values
   equivalent to epochal time epoch    */
void
etoh(Time *t, double epoch)
{
    int             diy;
    double          secleft;
    int            days;
    void            month_day();

    days = (int) (epoch / 86400.);
    secleft = mod(epoch, 86400.0);
    t->hr = t->mn = 0;
    t->sec = 0.;

    if (secleft) {		/* compute hours minutes seconds */
	if (secleft < 0) {	/* before 1970 */
	    days--;		/* subtract a day */
	    secleft += 86400.;	/* add a day */
	}
	t->hr = (int) (secleft / 3600.);
	secleft = mod(secleft, 3600.0);
	t->mn = (int) (secleft / 60.);
	t->sec = mod(secleft, 60.0);
    }

    if (days >= 0) {
	for (t->yr = 1970;; t->yr++) {
	    diy = isleap(t->yr) ? 366 : 365;
	    if (days < diy)
		break;
	    days -= diy;
	}
    }
    else {
	for (t->yr = 1969;; t->yr--) {
	    diy = isleap(t->yr) ? 366 : 365;
	    days += diy;
	    if (days >= 0)
		break;
	}
    }
    days++;
    month_day(t, (int) days);
    return;
}

/* set month and day fields of time structure *t
   given the julian day jul_day               */
void
month_day(Time *t, int jul_day)
{
    int             i, dim, leap;

    leap = isleap(t->yr);
    t->day = jul_day;
    for (i = 0; i < 12; i++) {
	dim = days_in_month[i];
	if (leap && (i == 1))
	    dim++;
	if (t->day <= dim)
	    break;
	t->day -= dim;
    }
    t->mo = i + 1;
    return;
}

int julian(Time *t)
/* returns the julian day represented by *t */
{
    int             i, j, inc;

    j = 0;
    for (i = 0; i < t->mo - 1; i++) {
	inc = days_in_month[i];
	if ((i == 1) && isleap(t->yr))
	    inc++;
	j += inc;
    }
    j += t->day;
    return (j);
}

int htoe(Time *t)
/* calculate epochal equivalent of time t */
{
    int            i, days;
    int            epoch;

    days = 0;
    if (t->yr > 1970) {
	for (i = 1970; i < t->yr; i++) {
	    days += 365;
	    if (isleap(i))
		days++;
	}
    }
    else if (t->yr < 1970) {
	for (i = t->yr; i < 1970; i++) {
	    days -= 365;
	    if (isleap(i))
		days--;
	}
    }
    days += (int) (julian(t) - 1);
    epoch = (((int) days) * 86400.);
    epoch += (int) ((t->hr * 3600.) + (t->mn * 60.) + (t->sec));
    return (epoch);
}

bool isleap(int year)
/* returns true if leap year, false otherwise */
{
    return ((year % 4 == 0 && year % 100 != 0) || year % 400 == 0);
}
void
output_(	float   data[],
		int    *ixs,
		float  *pas,
		int    *nt,
	        int    *dummy,
		char*   chan)
{
	char	file[132];
	int	i;
	Time	t;
	FILE    *fp;
	SAC     _sachead;
	SAC     *sachead;

	sachead=&_sachead;
	memcpy(sachead,&sac_null,sizeof(SAC));

/* Fake date */
        fp=fopen("axi.date","r");
        if (fp!=NULL) {
          fscanf(fp,"%d,%d,%d,%d,%d,%f",&t.yr,&t.mo,&t.day,&t.hr,&t.mn,&t.sec);
          fclose(fp);
	} 
        else
        {
          t.yr=2000;
          t.mo=1;
          t.day=1;
          t.hr=10;
          t.mn=0;
          t.sec=0.;
	}

	memset(sachead,0,sizeof(SAC));

	sachead->delta=*pas;
	sachead->b = 0;
	sachead->e = 0;
        sachead->nzyear = t.yr;
        sachead->nzjday = julian(&t);
        sachead->nzhour = t.hr;
        sachead->nzmin  = t.mn;
        sachead->nzsec  = (int)(t.sec);
        sachead->nzmsec = (t.sec-sachead->nzsec)*1000;

        sachead->nvhdr = 6;

        sachead->internal5 = 0;
        sachead->internal6 = 0;
        sachead->leven = TRUE;
        sachead->lpspol = FALSE;
        sachead->lcalda = TRUE;
        sachead->unused27 = FALSE;
        sachead->iftype = ITIME;
        sachead->iztype = IB;

	sprintf(sachead->kstnm,"%03d",*ixs);
	sprintf(sachead->kcmpnm,"%c",chan[0]);
	if (chan[0]=='X') {
		sachead->cmpinc=90.;
		sachead->cmpaz=0.;
	} else if (chan[0]=='Y'){
		sachead->cmpinc=90.;
		sachead->cmpaz=90.;
	} else if (chan[0]=='Z'){
		sachead->cmpinc=0.;
		sachead->cmpaz=0.;
	}

	sachead->npts = *nt;
	sachead->e=sachead->npts*sachead->delta;

	sprintf(file,"axi%03d.%c.sac",*ixs,chan[0]);
	fp=fopen(file,"w");
	fwrite(sachead,sizeof(SAC),1,fp);
	fwrite(data,sizeof(float),*nt,fp);
	fclose(fp);
}
