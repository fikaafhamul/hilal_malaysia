import streamlit as st
import math
from math import sin, cos, asin, degrees, radians, atan2, log10 , exp, sqrt

from pytz import timezone, utc
from datetime import timedelta

from skyfield.api import load, Topos
from skyfield.almanac import find_discrete, sunrise_sunset, moon_phases
from skyfield.nutationlib import iau2000b
from skyfield.framelib import ecliptic_frame

import pandas as pd

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.feature as cfeature
from cartopy import crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import seaborn as sns

import warnings
warnings.filterwarnings("ignore")

ts = load.timescale()
e = load('de440s.bsp')  # Menggunakan ephemeris DE440s
sun, moon, earth = e['sun'], e['moon'], e['earth']

class konversi:
    def __init__(self, angle, pilihan="DERAJAT"):
        self.angle = angle  # Menyimpan nilai angle ke dalam atribut instance
        pilihan = pilihan.upper()
        if pilihan == "DERAJAT":
            self.result = self.deg2dms()
        elif pilihan == "JAM":
            self.result = self.deg2hms()
        elif pilihan == "DERAJAT1":
            self.result = self.deg2dms1()
        elif pilihan == "JAM1":
            self.result = self.deg2hms1()
        elif pilihan == "LINTANG":
            self.result = self.deg2dms2()
        elif pilihan == "BUJUR":
            self.result = self.deg2dms3()
        else:
            self.result = self.deg2dms()

    def deg2dms(self):
        D = int(abs(self.angle))
        M = int((abs(self.angle) - D) * 60)
        S = (abs(self.angle) - D - M / 60) * 3600

        # Check for rounding issues
        if S >= 60:
            S = 0
            M += 1
        if M >= 60:
            M = 0
            D += 1

        if self.angle < 0:
            return f"-{D}° {M:02}' {S:05.2f}\""
        else:
            return f"{D}° {M:02}' {S:05.2f}\""

    def deg2hms(self):
        D = int(abs(self.angle))
        M = int((abs(self.angle) - D) * 60)
        S = (abs(self.angle) - D - M / 60) * 3600

        # Check for rounding issues
        if S >= 60:
            S = 0
            M += 1
        if M >= 60:
            M = 0
            D += 1

        if self.angle < 0:
            return f"-{D:02}:{M:02}:{S:05.2f}"
        else:
            return f"{D:02}:{M:02}:{S:05.2f}"

    def deg2dms1(self):
        D = int(abs(self.angle))
        M = int((abs(self.angle) - D) * 60)
        S = (abs(self.angle) - D - M / 60) * 3600

        # Check for rounding issues
        if S >= 60:
            S = 0
            M += 1
        if M >= 60:
            M = 0
            D += 1

        if self.angle < 0:
            if D == 0 and M == 0:
                return f"-{S:05.2f}\""
            elif D == 0:
                return f"-{M:02}' {S:05.2f}\""
            else:
                return f"-{D}° {M:02}' {S:05.2f}\""
        else:
            if D == 0 and M == 0:
                return f"{S:05.2f}\""
            elif D == 0:
                return f"{M:02}' {S:05.2f}\""
            else:
                return f"{D}° {M:02}' {S:05.2f}\""

    def deg2hms1(self):
        D = int(abs(self.angle))
        M = int((abs(self.angle) - D) * 60)
        S = (abs(self.angle) - D - M / 60) * 3600

        # Check for rounding issues
        if S >= 60:
            S = 0
            M += 1
        if M >= 60:
            M = 0
            D += 1

        if self.angle < 0:
            if D == 0 and M == 0:
                return f"-{S:05.2f} detik"
            elif D == 0:
                return f"-{M:02} Menit {S:05.2f} Detik"
            else:
                return f"-{D } jam {M:02} Menit {S:05.2f} Detik"
        else:
            if D == 0 and M == 0:
                return f"{S:05.2f} detik"
            elif D == 0:
                return f"{M:02} Menit {S:05.2f} Detik"
            else:
                return f"{D} jam {M:02} Menit {S:05.2f} Detik"
    def deg2dms2(self):
        D = int(abs(self.angle))
        M = int((abs(self.angle) - D) * 60)
        S = (abs(self.angle) - D - M / 60) * 3600

        # Check for rounding issues
        if S >= 60:
            S = 0
            M += 1
        if M >= 60:
            M = 0
            D += 1

        if self.angle < 0:
            return f"{D}° {M:02}' {S:05.2f}\" LS"
        else:
            return f"{D}° {M:02}' {S:05.2f}\" LU"

    def deg2dms3(self):
        D = int(abs(self.angle))
        M = int((abs(self.angle) - D) * 60)
        S = (abs(self.angle) - D - M / 60) * 3600

        # Check for rounding issues
        if S >= 60:
            S = 0
            M += 1
        if M >= 60:
            M = 0
            D += 1

        if self.angle < 0:
            return f"{D}° {M:02}' {S:05.2f}\" BB"
        else:
            return f"{D}° {M:02}' {S:05.2f}\" BT"

class hijriah:
    def __init__(self):
        self.bulan_hijriah = lambda bulan: self.b2_hijri(bulan) if isinstance(bulan, int) else self.hijri_to_b(bulan)

    def b2_hijri(self, bulan: int) -> str:
        if bulan == 1:
            return "Muharram"
        elif bulan == 2:
            return "Safar"
        elif bulan == 3:
            return "Rabiulawal"
        elif bulan == 4:
            return "Rabiulakhir"
        elif bulan == 5:
            return "Jamadilawal"
        elif bulan == 6:
            return "Jamadilakhir"
        elif bulan == 7:
            return "Rejab"
        elif bulan == 8:
            return "Syaaban"
        elif bulan == 9:
            return "Ramadan"
        elif bulan == 10:
            return "Syawal"
        elif bulan == 11:
            return "Zulkaedah"
        elif bulan == 12:
            return "Zulhijjah"
        else:
            return ""

    def hijri_to_b(self, bulan: str) -> int:
        if bulan == "Muharram":
            return 1
        elif bulan == "Safar":
            return 2
        elif bulan == "Rabiulawal":
            return 3
        elif bulan == "Rabiulakhir":
            return 4
        elif bulan == "Jumadilawal":
            return 5
        elif bulan == "Jumadilakhir":
            return 6
        elif bulan == "Rejab":
            return 7
        elif bulan == "Syaaban":
            return 8
        elif bulan == "Ramadan":
            return 9
        elif bulan == "Syawal":
            return 10
        elif bulan == "Zulkaedah":
            return 11
        elif bulan == "Zulhijjah":
            return 12
        else:
            return 0

class miladi:
    def __init__(self):
        self.bulan_miladi = lambda bulan: self.b2_miladi(bulan) if isinstance(bulan, int) else self.miladi_to_b(bulan)

    def b2_miladi(self, bulan: int) -> str:
        if bulan == 1:
            return "Januari"
        elif bulan == 2:
            return "Februari"
        elif bulan == 3:
            return "Mac"
        elif bulan == 4:
            return "April"
        elif bulan == 5:
            return "Mei"
        elif bulan == 6:
            return "Jun"
        elif bulan == 7:
            return "Julai"
        elif bulan == 8:
            return "Ogos"
        elif bulan == 9:
            return "September"
        elif bulan == 10:
            return "Oktober"
        elif bulan == 11:
            return "November"
        elif bulan == 12:
            return "Disember"
        else:
            return ""

    def miladi_to_b(self, bulan: str) -> int:
        if bulan == "Januari":
            return 1
        elif bulan == "Februari":
            return 2
        elif bulan == "Mac":
            return 3
        elif bulan == "April":
            return 4
        elif bulan == "Mei":
            return 5
        elif bulan == "Jun":
            return 6
        elif bulan == "Julai":
            return 7
        elif bulan == "Ogos":
            return 8
        elif bulan == "September":
            return 9
        elif bulan == "Oktober":
            return 10
        elif bulan == "November":
            return 11
        elif bulan == "Disember":
            return 12
        else:
            return 0

class caldat:

    def __init__(self, jd, timezone=0.0, pilihan=None):
        self.jd = jd
        self.timezone = timezone
        self.pilihan = pilihan
        self.hari_str, self.pasaran_str, self.tgl, self.bln, self.thn, self.jam = self.calculate()

        # Memastikan pilihan tidak None dan mengubahnya menjadi uppercase
        pilihan = self.pilihan.upper() if self.pilihan else None

        # Mengatur hasil berdasarkan pilihan
        if pilihan == "HARI":
            self.result = self.hari_str
        elif pilihan == "PASARAN":
            self.result = self.pasaran_str
        elif pilihan == "HARPAS":
            self.result = f"{self.hari_str} {self.pasaran_str}"
        elif pilihan == "TANGGAL":
            self.result = f"{self.tgl:02d} {miladi().bulan_miladi(self.bln)} {self.thn:04d}"
        elif pilihan == "JAM":
            self.result = konversi(self.jam, "JAM").result()
        elif pilihan == "JDJAM":
            self.result = self.jam
        elif pilihan == "JDTANGGAL":
            self.result = self.tgl
        elif pilihan == "JDBULAN":
            self.result = self.bln
        elif pilihan == "JDTAHUN":
            self.result = self.thn
        elif pilihan == "JD_LENGKAP":
            self.result = self.tgl, self.bln, self.thn, self.hari_str, self.pasaran_str
        elif pilihan == "JD_HP":
            self.result = f"{self.hari_str}\t{self.pasaran_str}"
        else:
            self.result = f"{self.hari_str}, {self.tgl:02d} {miladi().bulan_miladi(self.bln)} {self.thn:04d}"


    def calculate(self):
        z = int(self.jd + 0.5)

        if z < 2299161:
            a = z
        else:
            alpha = int((z - 1867216.25) / 36524.25)
            a = z + 1 + alpha - int(alpha / 4)

        jam = ((self.jd + 0.5 + self.timezone / 24) - z) * 24
        if jam > 24:
            jam -= 24
            a += 1
            z += 1

        b = a + 1524
        c = int((b - 122.1) / 365.25)
        d = int(365.25 * c)
        e = int((b - d) / 30.6001)

        tgl = int(b - d - int(30.6001 * e))

        if e < 14:
            bln = e - 1
        elif e == 14 or e == 15:
            bln = e - 13

        if bln > 2:
            thn = c - 4716
        else:
            thn = c - 4715

        hari = math.fmod(int(z) + 2, 7)
        hari_nama = ["Sabtu", "Ahad", "Isnin", "Selasa", "Rabu", "Khamis", "Jumaat"]
        hari_str = hari_nama[int(hari)]

        pasaran = math.fmod(int(z) + 1, 5)
        pasaran_nama = ["Kliwon", "Legi", "Pahing", "Pon", "Wage"]
        pasaran_str = pasaran_nama[int(pasaran)]

        return hari_str, pasaran_str, tgl, bln, thn, jam

class awalbulan:
    def __init__(self, bulan, tahun, lok, lat, lon, TZ='Asia/Kuala_Lumpur', TT = 0, TH = 0, kriteria = 'NEO MABIMS'):
        self.bulan = bulan
        self.tahun = tahun
        if self.bulan < 2:
            self.bulan1 = bulan - 1 + 12
            self.tahun1 = tahun - 1
        else:
            self.bulan1 = bulan - 1
            self.tahun1 = tahun
        self.lokasi = lok  # Nama Lokasi
        self.lat = lat  # Latitude pengamat
        self.lon = lon  # Longitude pengamat
        self.TZ = TZ
        self.TT = TT
        self.TH = TH
        self.kriteria = kriteria
        self.JDE = self.hitung_jde()  # Menghitung JDE saat inisialisasi
        self.newmoon = self.new_moon()  # Mengambil nilai konjungsi saat inisialisasi
        self.moonrise_moonset = self.rise_set_moon()

    def hitung_jde(self):
        # Menghitung Hy
        Hy = self.tahun1 + (((self.bulan1 - 1) * 29.53) / 354.3671)

        # Menghitung K
        K = round(((Hy - 1410) * 12), 0) - 129

        # Menghitung T
        T = K / 1236.85

        # Menghitung JDE
        JDE = 2451550.09765 + 29.530588853 * K + 0.0001337 * (T ** 2)

        return JDE

    def new_moon(self):
        temp = caldat(self.JDE, 0, "JD_LENGKAP").result
        temp = list(temp)

        # Mengonversi JDE ke waktu UTC
        jd_time = ts.utc(temp[2], temp[1], temp[0], 0)

        # Mengatur rentang waktu untuk pencarian New Moon
        t0 = jd_time - timedelta(days=2)
        t1 = jd_time + timedelta(days=2)

        # Mencari fase bulan
        phases_time, y = find_discrete(t0, t1, moon_phases(e))

        # Mengambil waktu New Moon
        new_moon_times = phases_time[y == 0]  # 0 menandakan New Moon

        # Mengatur zona waktu (WIB)
        ZonaWaktu = timezone(self.TZ)
        return new_moon_times.astimezone(ZonaWaktu)

    def rise_set_moon(self):
        moon = e['moon']
        longlat = Topos(latitude=self.lat, longitude=self.lon)
        topos_at = (e['earth'] + longlat).at
        def is_moon_up_at(t):
            t._nutation_angles = iau2000b(t.tt)
            return topos_at(t).observe(moon).apparent().altaz()[0].degrees > -50/60
        is_moon_up_at.rough_period = 0.5
        return is_moon_up_at

    def calculate_hilal(self):
        # Ambil waktu konjungsi pertama
        konjungsi_time = self.newmoon[0]

        # Menentukan rentang waktu untuk pencarian matahari terbenam
        t0 = ts.utc(konjungsi_time.year, konjungsi_time.month, konjungsi_time.day, 0, 0)
        t1 = t0 + timedelta(days=1)

        # Menghitung Lintang dan Bujur Pengamat
        longlat = Topos(latitude=self.lat, longitude=self.lon, elevation_m=self.TT)

        # Menghitung waktu terbenam matahari
        sunriset, sunBol = find_discrete(t0, t1, sunrise_sunset(e, longlat))
        sunset_time = sunriset[sunBol == 0]  # 0 menandakan waktu terbenam

        # Ubah ke waktu UTC
        sunset_time_utc = sunset_time.utc_iso()

        # Ubah ke waktu lokal
        ZonaWaktu = timezone(self.TZ)
        sunset_time_local = sunset_time.astimezone(ZonaWaktu)

        # Menghitung Umur Bulan
        temp = konjungsi_time.hour + (konjungsi_time.minute)/60 + (konjungsi_time.second)/3600
        temp1 = sunset_time_local[0]
        temp1 = temp1.hour + (temp1.minute)/60 + (temp1.second)/3600
        moonage = (temp1 - temp)

        # Menghitung waktu terbenam Bulan
        moonriset, moonBol = find_discrete(t0, t1, self.moonrise_moonset)
        moonset_time = moonriset[moonBol == 0]

        # Ubah ke waktu UTC
        moonset_time_utc = moonset_time.utc_iso()

        # Ubah ke waktu lokal
        moonset_time_local = moonset_time.astimezone(ZonaWaktu)

        observer = earth + longlat

        # Menghitung Tinggi dan Elongasi Bulan
        geo_moon = earth.at(sunset_time[0]).observe(moon).apparent()
        geo_sun = earth.at(sunset_time[0]).observe(sun).apparent()
        el_geo = geo_sun.separation_from(geo_moon).degrees

        topo_moon = observer.at(sunset_time[0]).observe(moon).apparent()
        topo_sun = observer.at(sunset_time[0]).observe(sun).apparent()
        alt, az, distance = topo_moon.altaz()
        alt = alt.degrees
        el_topo = topo_sun.separation_from(topo_moon).degrees

        kriteria = self.kriteria.upper()
        if kriteria == "IRNU":
            jd = (t0 + timedelta(days=1)) if (alt >= 3 and el_topo >= 6.4) or (el_topo > 9.9) else (t0 + timedelta(days=2))
        elif kriteria == "MUHAMADIYYAH":
            jd = (t0 + timedelta(days=1)) if (alt >= 0) else (t0 + timedelta(days=2))
        elif kriteria == "MABIMS LAMA":
            jd = (t0 + timedelta(days=1)) if (alt >= 2 and el_topo >= 3 and moonage >= 8) else (t0 + timedelta(days=2))
        elif kriteria == "NEO MABIMS":
            jd = (t0 + timedelta(days=1)) if (alt >= 3 and el_topo >= 6.4) else (t0 + timedelta(days=2))

        jd = jd + timedelta(days=28) + self.TH

        # Mengatur rentang waktu untuk pencarian New Moon
        t0 = jd - timedelta(days=2)
        t1 = jd + timedelta(days=2)

        # Mencari fase bulan
        phases_time, y = find_discrete(t0, t1, moon_phases(e))

        # Mengambil waktu New Moon
        new_moon_times = phases_time[y == 0]  # 0 menandakan New Moon
        new_moon_times = new_moon_times.astimezone(ZonaWaktu)

        konjungsi_times = new_moon_times[0]

        sunriset, sunBol = find_discrete(jd, jd + timedelta(days=1), sunrise_sunset(e, longlat))
        sunset_time = sunriset[sunBol == 0]  # 0 menandakan waktu terbenam

        # Ubah ke waktu UTC
        sunset_time_utc = sunset_time.utc_iso()

        # Ubah ke waktu lokal
        sunset_time_local = sunset_time.astimezone(ZonaWaktu)

        # Menghitung Umur Bulan
        temp = konjungsi_times.hour + (konjungsi_times.minute)/60 + (konjungsi_times.second)/3600
        temp1 = sunset_time_local[0]
        temp1 = temp1.hour + (temp1.minute)/60 + (temp1.second)/3600
        moonage = (temp1 - temp)

        # Menghitung waktu terbenam Bulan
        moonriset, moonBol = find_discrete(jd, jd + timedelta(days=1), self.moonrise_moonset)
        moonset_time = moonriset[moonBol == 0]

        # Ubah ke waktu UTC
        moonset_time_utc = moonset_time.utc_iso()

        # Ubah ke waktu lokal
        moonset_time_local = moonset_time.astimezone(ZonaWaktu)

        # Menghitung Tinggi dan Elongasi Bulan
        geo_moon = earth.at(sunset_time[0]).observe(moon).apparent()
        geo_sun = earth.at(sunset_time[0]).observe(sun).apparent()
        el_geo = geo_sun.separation_from(geo_moon).degrees

        topo_moon = observer.at(sunset_time[0]).observe(moon).apparent()
        topo_sun = observer.at(sunset_time[0]).observe(sun).apparent()
        moon_alt, moon_az, moon_distance = topo_moon.altaz()
        moon_alt, moon_az = moon_alt.degrees, moon_az.degrees
        sun_alt, sun_az, sun_distance = topo_sun.altaz()
        sun_alt, sun_az = sun_alt.degrees, sun_az.degrees
        el_topo = topo_sun.separation_from(topo_moon).degrees

        # Menghitung Iluminasi Bulan
        illumination = topo_moon.fraction_illuminated(sun)*100

        # Menghitung Lebar Sabit Bulan
        horizontalparalax = degrees(asin(6378.14/moon_distance.km))
        semidiameter= (358473400/moon_distance.km)/3600
        SD = semidiameter * horizontalparalax
        SD1 = SD * (1 + sin(radians(moon_alt))*sin(radians(horizontalparalax)))
        W = SD1 * (1-cos(radians(el_topo)))

        # Menghitung Umur Bulan
        temp = konjungsi_times.hour + (konjungsi_times.minute)/60 + (konjungsi_times.second)/3600
        temp1 = sunset_time_local[0]
        temp1 = temp1.hour + (temp1.minute)/60 + (temp1.second)/3600
        moonage = (temp1 - temp)

        sunset = sunset_time_local[0]
        moonset = moonset_time_local[0]

        # Lama Hilal
        temp = sunset.hour + (sunset.minute)/60 + (sunset.second)/3600
        temp1 = moonset.hour + (moonset.minute)/60 + (moonset.second)/3600
        lag_time = (temp1 - temp)

        separasi = moon_alt - sun_alt
        best_time, q, parameter = visibilitas_oddeh(temp, lag_time, separasi, W)

        n_bln = miladi().bulan_miladi(sunset.month)
        temp = sunset.tzinfo.utcoffset(sunset)
        delta_time_tz = int(temp.total_seconds()/3600)

        return {
            "konjungsi": konjungsi_times,
            "jd": jd,
            "sunset": sunset,
            "moonset": moonset,
            "sun_alt": sun_alt,
            "sun_az": sun_az,
            "moon_alt": moon_alt,
            "moon_az": moon_az,
            "moonage": moonage,
            "el_geo": el_geo,
            "el_topo": el_topo,
            "moon_age": moonage,
            "illumination": illumination,
            "lag_time": lag_time,
            "lebar_sabit": W,
            "best_time": best_time,
            "q": q,
            "parameter": parameter,
            "dt": delta_time_tz
        }

def visibilitas_oddeh(sunset, lag, separasi, lebar_sabit):
    # Menghitung Best Time
    best_time = sunset + 4/9 * lag

    # Menghitung q
    lebar_sabit = lebar_sabit * 60
    q = separasi - (-0.1018 * (lebar_sabit ** 3) + 0.7319 * (lebar_sabit ** 2) - 6.3226 * lebar_sabit + 7.1651)

    # Parameter
    if q >= 5.65:
        parameter = "Mudah diperhatikan tanpa bantuan optik"
    elif (q < 5.65 and q >= 2):
        parameter = "Boleh diperhatikan tanpa bantuan optik"
    elif (q < 2 and q >= -0.96):
        parameter = "Kelihatan hanya dengan bantuan optik"
    elif (q < -0.96):
        parameter = "Tak nampak"

    return best_time, q, parameter

def visibilitas_katsner(AwalBulan, latitude, longitude, TZ, TT):
    df = pd.DataFrame(columns=['T AFTER SUNSET', 'TIME', 'SUN ALTITUDE', 'MOON ALTITUDE', 'SUN AZIMUT', 'MOON AZIMUT','MOON ELONGATION',
                               'NAKED EYE (K2)', 'NAKED EYE (K4)', 'NAKED EYE (K6)', 'THEODOLIT (K2)', 'THEODOLIT (K4)', 'THEODOLIT (K6)',
                               'TELESCOPE (K2)', 'TELESCOPE (K4)', 'TELESCOPE (K6)'])
    ZonaWaktu = timezone(TZ)

    data = variabel(AwalBulan)

    longlat = Topos(latitude=latitude, longitude=longitude, elevation_m=TT)
    observer = earth + longlat

    sunset = data["sunset"]
    moonset = data["moonset"]
    sunset = sunset.astimezone(utc)
    moonset = moonset.astimezone(utc)
    t = ts.utc(sunset.year, sunset.month, sunset.day, sunset.hour, sunset.minute)
    t_end = t + timedelta(hours=1)

    t_moonset = ts.utc(moonset.year, moonset.month, moonset.day, moonset.hour, moonset.minute)
    T_after_sunset = 0

    def koreksi_theodolit():
        #Fb
        Fb = sqrt(2)

        #Fr
        M = 30
        theta_s = 3
        M_Theta = 2 * M *theta_s
        Fr = (M_Theta / 900) ** 0.5 if M_Theta > 900 else 1

        #Fp
        usia_pengamat = 23
        diameter_pupil = 7 * exp(-0.5 * (usia_pengamat / 100) ** 2)
        diameter_objektif = 45
        diameter_exit_pupil = diameter_objektif / M
        Fp = (diameter_objektif / (M * diameter_pupil)) ** 2 if diameter_pupil < diameter_exit_pupil else 1

        #Ft
        sigma_permukaan = 6
        faktor_transmisi = 0.95
        Ds = diameter_pupil / 71
        Ft = 1/((faktor_transmisi ** sigma_permukaan) * (1 - (Ds) ** 2))

        #Fa
        Fa = (diameter_pupil / diameter_objektif) ** 2

        #Fm
        Fm = M ** 2

        return Fb, Fr, Fp, Ft, Fa, Fm

    def koreksi_telescope():
        #Fb
        Fb = sqrt(2)

        #Fr
        M = 32
        theta_s = 3
        M_Theta = 2 * M *theta_s
        Fr = (M_Theta / 900) ** 0.5 if M_Theta > 900 else 1

        #Fp
        usia_pengamat = 23
        diameter_pupil = 7 * exp(-0.5 * (usia_pengamat / 100) ** 2)
        diameter_objektif = 80
        diameter_exit_pupil = diameter_objektif / M
        Fp = (diameter_objektif / (M * diameter_pupil)) ** 2 if diameter_pupil < diameter_exit_pupil else 1

        #Ft
        sigma_permukaan = 6
        faktor_transmisi = 0.95
        Ds = diameter_pupil / 71
        Ft = 1/((faktor_transmisi ** sigma_permukaan) * (1 - (Ds) ** 2))

        #Fa
        Fa = (diameter_pupil / diameter_objektif) ** 2

        #Fm
        Fm = M ** 2

        return Fb, Fr, Fp, Ft, Fa, Fm

    while t.utc_datetime() <= t_end.utc_datetime():

        geo_moon = earth.at(t).observe(moon).apparent()
        geo_sun = earth.at(t).observe(sun).apparent()
        el_geo = geo_sun.separation_from(geo_moon).degrees

        topo_moon = observer.at(t).observe(moon).apparent()
        topo_sun = observer.at(t).observe(sun).apparent()
        moon_alt, moon_az, moon_distance = topo_moon.altaz()
        moon_alt, moon_az = moon_alt.degrees, moon_az.degrees
        sun_alt, sun_az, sun_distance = topo_sun.altaz()
        sun_alt, sun_az = sun_alt.degrees, sun_az.degrees
        el_topo = topo_sun.separation_from(topo_moon).degrees

        if t.utc_datetime() > t_moonset.utc_datetime():
            break

        moon_zenit = 90 - moon_alt
        sun_depression = abs(sun_alt)
        rel_azimut = abs(moon_az - sun_az)
        O_zero = -(4.12 * 10 ** (-2) * moon_zenit + 0.582) * sun_depression + 0.417 * moon_zenit + 97.5
        log_L = -(7.5 * 10 ** (-5) * moon_zenit + 5.05 * 10 ** (-3)) * rel_azimut + (3.67 * 10 ** (-4) * moon_zenit - 0.458) * sun_depression + 9.17 * 10 ** (-3) * moon_zenit + 3.525
        log_L1 = -(7.5 * 10 ** (-5) * 90 + 5.05 * 10 ** (-3)) * rel_azimut + (3.67 * 10 ** (-4) * 90 - 0.458) * sun_depression + 9.17 * 10 ** (-3) * 90 + 3.525
        log_L2 = -(7.5 * 10 ** (-5) * 30 + 5.05 * 10 ** (-3)) * rel_azimut + (3.67 * 10 ** (-4) * 30 - 0.458) * sun_depression + 9.17 * 10 ** (-3) * 30 + 3.525
        log_LI = (((moon_zenit - 60) / 30) * log_L1) + (((90 - moon_zenit) / 30) * log_L2)
        moon_mvis = 0.026 * degrees(atan2(sun_distance.km * sin(radians(el_geo)) , moon_distance.km - sun_distance.km * cos(radians(el_geo)))) + ((4 * 10 ** (-9)) * (degrees(atan2(sun_distance.km * sin(radians(el_geo)) , moon_distance.km - sun_distance.km * cos(radians(el_geo)))) ** 4)) - 12.73
        moon_semidiameter= (358473400/moon_distance.km)/3600
        crescent_square_degrees = 0.5 * (22/7) * (moon_semidiameter ** 2) * (1 + cos(radians(180-el_topo)))
        crescent_square_arc = crescent_square_degrees * 3600 * 3600
        moon_SB = moon_mvis + (2.5 * log10(crescent_square_arc))
        moon_L_extra = (1 / crescent_square_degrees) * 2.51 ** (10 - moon_mvis)
        term_1st = cos(radians(moon_zenit))
        term_2nd = 0.025 * exp(-11 * cos(radians(moon_zenit)))
        x_final = 1 / (term_1st + term_2nd)

        ls_S10_calibrated = 290 * (10 ** (log_LI + 2.5))
        la_S10 = 290 + 105 * exp(-(90 - moon_zenit) ** 2 / 1600)

        def naked_eye():

            # Atmosfer dengan k = 0,2
            moon_L_ground = moon_L_extra * exp(-0.2 * x_final)
            R = moon_L_ground / (ls_S10_calibrated + la_S10)
            Delta_M_K2 = 2.5 * log10(R)

            # Atmosfer dengan k = 0,4
            moon_L_ground = moon_L_extra * exp(-0.4 * x_final)
            R = moon_L_ground / (ls_S10_calibrated + la_S10)
            Delta_M_K4 = 2.5 * log10(R)

            # Atmosfer dengan k = 0,6
            moon_L_ground = moon_L_extra * exp(-0.6 * x_final)
            R = moon_L_ground / (ls_S10_calibrated + la_S10)
            Delta_M_K6 = 2.5 * log10(R)
            return {"K2": Delta_M_K2,
                    "K4": Delta_M_K4,
                    "K6": Delta_M_K6}

        def theodolite():
            Fb, Fr, Fp, Ft, Fa, Fm = koreksi_theodolit()

            # Atmosfer dengan k = 0,2
            moon_L_ground = moon_L_extra * exp(-0.2 * x_final)
            moon_L_ground_teleskop = moon_L_ground / (Fb * Fr * Fp * Ft * Fa)
            log_LI_teleskop = log10((ls_S10_calibrated + la_S10) / (Fb * Fp * Ft * Fa * Fm))
            ls_S10_calibrated_teleskop = 10 ** log_LI_teleskop
            R = moon_L_ground_teleskop / ls_S10_calibrated_teleskop
            Delta_M_K2 = 2.5 * log10(R)

            # Atmosfer dengan k = 0,4
            moon_L_ground = moon_L_extra * exp(-0.4 * x_final)
            moon_L_ground_teleskop = moon_L_ground / (Fb * Fr * Fp * Ft * Fa)
            log_LI_teleskop = log10((ls_S10_calibrated + la_S10) / (Fb * Fp * Ft * Fa * Fm))
            ls_S10_calibrated_teleskop = 10 ** log_LI_teleskop
            R = moon_L_ground_teleskop / ls_S10_calibrated_teleskop
            Delta_M_K4 = 2.5 * log10(R)

            # Atmosfer dengan k = 0,6
            moon_L_ground = moon_L_extra * exp(-0.6 * x_final)
            moon_L_ground_teleskop = moon_L_ground / (Fb * Fr * Fp * Ft * Fa)
            log_LI_teleskop = log10((ls_S10_calibrated + la_S10) / (Fb * Fp * Ft * Fa * Fm))
            ls_S10_calibrated_teleskop = 10 ** log_LI_teleskop
            R = moon_L_ground_teleskop / ls_S10_calibrated_teleskop
            Delta_M_K6 = 2.5 * log10(R)

            return {"K2": Delta_M_K2,
                    "K4": Delta_M_K4,
                    "K6": Delta_M_K6}

        def telescope():
            Fb, Fr, Fp, Ft, Fa, Fm = koreksi_telescope()

            # Atmosfer dengan k = 0,2
            moon_L_ground = moon_L_extra * exp(-0.2 * x_final)
            moon_L_ground_teleskop = moon_L_ground / (Fb * Fr * Fp * Ft * Fa)
            log_LI_teleskop = log10((ls_S10_calibrated + la_S10) / (Fb * Fp * Ft * Fa * Fm))
            ls_S10_calibrated_teleskop = 10 ** log_LI_teleskop
            R = moon_L_ground_teleskop / ls_S10_calibrated_teleskop
            Delta_M_K2 = 2.5 * log10(R)

            # Atmosfer dengan k = 0,4
            moon_L_ground = moon_L_extra * exp(-0.4 * x_final)
            moon_L_ground_teleskop = moon_L_ground / (Fb * Fr * Fp * Ft * Fa)
            log_LI_teleskop = log10((ls_S10_calibrated + la_S10) / (Fb * Fp * Ft * Fa * Fm))
            ls_S10_calibrated_teleskop = 10 ** log_LI_teleskop
            R = moon_L_ground_teleskop / ls_S10_calibrated_teleskop
            Delta_M_K4 = 2.5 * log10(R)

            # Atmosfer dengan k = 0,6
            moon_L_ground = moon_L_extra * exp(-0.6 * x_final)
            moon_L_ground_teleskop = moon_L_ground / (Fb * Fr * Fp * Ft * Fa)
            log_LI_teleskop = log10((ls_S10_calibrated + la_S10) / (Fb * Fp * Ft * Fa * Fm))
            ls_S10_calibrated_teleskop = 10 ** log_LI_teleskop
            R = moon_L_ground_teleskop / ls_S10_calibrated_teleskop
            Delta_M_K6 = 2.5 * log10(R)

            return {"K2": Delta_M_K2,
                    "K4": Delta_M_K4,
                    "K6": Delta_M_K6}

        naked_eye_delta_m = naked_eye()
        theodolit_delta_m = theodolite()
        telescope_delta_m = telescope()

        waktu = t.astimezone(ZonaWaktu)
        waktu = waktu.strftime("%H:%M")

        df_baru = pd.DataFrame({'T AFTER SUNSET' : [T_after_sunset],
                                'TIME' : [waktu],
                                'SUN ALTITUDE' : [konversi(sun_alt).result],
                                'MOON ALTITUDE' : [konversi(moon_alt).result],
                                'SUN AZIMUT' : [konversi(sun_az).result],
                                'MOON AZIMUT' : [konversi(moon_az).result],
                                'MOON ELONGATION' : [konversi(el_topo).result],
                                'NAKED EYE (K2)' : [naked_eye_delta_m["K2"]],
                                'NAKED EYE (K4)' : [naked_eye_delta_m["K4"]],
                                'NAKED EYE (K6)' : [naked_eye_delta_m["K6"]],
                                'THEODOLIT (K2)' : [theodolit_delta_m["K2"]],
                                'THEODOLIT (K4)' : [theodolit_delta_m["K4"]],
                                'THEODOLIT (K6)' : [theodolit_delta_m["K6"]],
                                'TELESCOPE (K2)' : [telescope_delta_m["K2"]],
                                'TELESCOPE (K4)' : [telescope_delta_m["K4"]],
                                'TELESCOPE (K6)' : [telescope_delta_m["K6"]]})

        t = t + timedelta(minutes=1)
        T_after_sunset += 1


        df = pd.concat([df, df_baru], ignore_index=True)

    return df

def cari_data(df):
    iterasi = ['NAKED EYE (K2)','NAKED EYE (K4)','NAKED EYE (K6)', 'THEODOLIT (K2)', 'THEODOLIT (K4)', 'THEODOLIT (K6)', 'TELESCOPE (K2)', 'TELESCOPE (K4)', 'TELESCOPE (K6)']
    first_index = []
    best_index = []
    for i in iterasi:
      cek = df[i].max()
      if cek > 0:
          best_index.append(df[i].idxmax())
      else:
          best_index.append('-')
      for index, temp in df[i].items():
          if temp > 0:
              first_index.append(index)
              break
          if index == df[i].index[-1]:
              first_index.append('-')
    return{
        'first_index' : first_index,
        'best_index' : best_index
    }

def plot_katsner_curve(df, AwalBulan, bulan_hijriah, tahun_hijriah, alat='TELESCOPE'):

    jd_map = AwalBulan["sunset"]
    jd_map = ts.utc(jd_map.year, jd_map.month, jd_map.day, 0)
    jd_caldat_result = caldat(float(jd_map.tt), 7).result

    iterasi = ['NAKED EYE (K2)','NAKED EYE (K4)','NAKED EYE (K6)', 'THEODOLIT (K2)', 'THEODOLIT (K4)', 'THEODOLIT (K6)', 'TELESCOPE (K2)', 'TELESCOPE (K4)', 'TELESCOPE (K6)']
    temp1 = cari_data(df)
    first = []
    best = []
    curve = []
    for i in iterasi:
        if i.startswith(alat):
            if temp1['first_index'][iterasi.index(i)] == '-':
                first.append('-')
            else:
                first.append(df['TIME'][temp1['first_index'][iterasi.index(i)]])
            if temp1['best_index'][iterasi.index(i)] == '-':
                best.append('-')
                curve.append('-')
            else:
                best.append(df['TIME'][temp1['best_index'][iterasi.index(i)]])
                curve.append(df['T AFTER SUNSET'][temp1['best_index'][iterasi.index(i)]])

    fig, ax = plt.subplots(1,1, figsize=(10,6))

    sns.lineplot(x='T AFTER SUNSET', y=f'{alat} (K2)', data=df, color='blue', linestyle='-', marker = 'o', label='Atmosfera dengan k = 0,2 (Bersih)', ax=ax)
    sns.lineplot(x='T AFTER SUNSET', y=f'{alat} (K4)', data=df, color='orange', linestyle='-', marker = 'o', label='Atmosfera dengan k = 0,4 (Moderat)', ax=ax)
    sns.lineplot(x='T AFTER SUNSET', y=f'{alat} (K6)', data=df, color='red', linestyle='-',  marker = 'o', label='Atmosfera dengan k = 0,6 (Kotor)', ax=ax)

    ax.set_xlabel('Minit (Selepas Matahari Terbenam)', fontsize=12)
    ax.set_ylabel('Kebolehnampakan', fontsize=12)

    bulan_hijriah_str = hijriah().bulan_hijriah(bulan_hijriah)

    plt.title(f'KELUK KEBOLEHNAMPAKAN MODEL KASTNER \n PENCERAPAN ANAK BULAN {bulan_hijriah_str.upper()} {tahun_hijriah} H \n {jd_caldat_result.upper()} \n', fontsize=14, pad = 20)

    ax.set_ylim(-30.5, 30)

    ax.spines['bottom'].set_position(('data', 0))
    plt.legend(loc='best')
    plt.grid(color='lightgray')
    ax.spines['left'].set_color('lightgray')
    ax.spines['top'].set_color('lightgray')
    ax.spines['right'].set_color('lightgray')

    ax.tick_params(axis='y', length=0)
    ax.tick_params(axis='x', length=0)

    ax.xaxis.labelpad = 170

    colors = ['blue', 'orange','red']
    atm = ['k = 0.2', 'k = 0.4', 'k = 0.6']
    for index, i in enumerate(curve):
        if i != '-':
            plt.axvline(x=i, color=colors[index], ls='--', lw=1, label = f'Best Time {atm[index]}')

    # Adjust legend position
    ax.legend(loc='upper left', bbox_to_anchor=(0, -0.1))
    # Adjust figtext positions
    plt.figtext(0.46, -0.08, f"First Visibility :\n     Keadaan Bersih    : {first[0]}\n     Keadaan Moderat : {first[1]}\n     Keadaan Kotor      : {first[2]}", fontsize=10)
    plt.figtext(0.68, -0.08, f"Best Time :\n     Keadaan Bersih    : {best[0]}\n     Keadaan Moderat : {best[1]}\n     Keadaan Kotor      : {best[2]}", fontsize=10)
    plt.figtext(0.45, 0.915, f"MODUS: {alat}", fontsize=10)

    return fig

def variabel(AwalBulan):
    n1_bln = miladi().bulan_miladi(AwalBulan["konjungsi"].month)
    jd = AwalBulan["jd"]
    konjungsi = AwalBulan["konjungsi"]

    delta_time = AwalBulan["dt"]
    best_time = AwalBulan["best_time"]
    q = AwalBulan["q"]
    parameter = AwalBulan["parameter"]

    sunset = AwalBulan["sunset"]
    n_bln = miladi().bulan_miladi(sunset.month)
    sun_alt = AwalBulan["sun_alt"]
    sun_az = AwalBulan["sun_az"]

    moonset = AwalBulan["moonset"]
    moon_alt = AwalBulan["moon_alt"]
    moon_az = AwalBulan["moon_az"]
    el_geo = AwalBulan["el_geo"]
    el_topo = AwalBulan["el_topo"]
    moonage = AwalBulan["moonage"]
    lag_time = AwalBulan["lag_time"]
    W = AwalBulan["lebar_sabit"]
    illumination = AwalBulan["illumination"]

    return {
        "n1_bln": n1_bln,
        "jd": jd,
        "konjungsi": konjungsi,
        "delta_time": delta_time,
        "best_time": best_time,
        "q": q,
        "parameter": parameter,
        "sunset": sunset,
        "n_bln": n_bln,
        "sun_alt": sun_alt,
        "sun_az": sun_az,
        "moonset": moonset,
        "moon_alt": moon_alt,
        "moon_az": moon_az,
        "el_geo": el_geo,
        "el_topo": el_topo,
        "moonage": moonage,
        "lag_time": lag_time,
        "W": W,
        "illumination": illumination
    }

def datapeta(AwalBulan):
    # Ambil data sunset dari dictionary AwalBulan
    sunset = AwalBulan["sunset"]

    # Definisikan rentang lintang dan bujur untuk peta
    latitudes = np.linspace(1, 7, 10)  # Indonesia's latitude range
    longitudes = np.linspace(100, 120, 10)  # Indonesia's longitude range

    # Inisialisasi array untuk menyimpan hasil
    altitudes = np.zeros((len(latitudes), len(longitudes)))
    elongations = np.zeros((len(latitudes), len(longitudes)))

    # Definisikan zona waktu berdasarkan bujur
    def determine_time_zone(lon):
        # if lon < 112.5:
        #     return timezone('Asia/Jakarta')  # WIB (UTC+7)
        # elif 112.5 <= lon < 127.5:
        #     return timezone('Asia/Makassar')  # WITA (UTC+8)
        # else:
        #     return timezone('Asia/Jayapura')  # WIT (UTC+9)
        return timezone('Asia/Kuala_Lumpur')

    # Definisikan rentang waktu pencarian dalam UTC
    dates = ts.utc(sunset.year, sunset.month, sunset.day, 0)
    t0 = dates
    t1 = dates + 1

    # Iterasi melalui lokasi-lokasi pada grid
    for i, lat in enumerate(latitudes):
        for j, lon in enumerate(longitudes):
            observer = Topos(latitude_degrees=lat, longitude_degrees=lon)
            observer_at = earth + observer
            observer_at_geocenter = earth

            # Cari waktu matahari terbenam dalam UTC untuk lokasi ini
            f = sunrise_sunset(e, observer)
            sunset_times, is_sunsets = find_discrete(t0, t1, f)
            sunset_times = sunset_times[is_sunsets == False]  # Ambil hanya waktu terbenam

            # Tentukan zona waktu lokal dan sesuaikan waktu terbenam
            local_tz = determine_time_zone(lon)

            # Untuk setiap waktu terbenam yang ditemukan (biasanya hanya satu dalam rentang 1 hari)
            for sunset_time_utc in sunset_times:
                 # Hitung ketinggian dan elongasi Bulan saat matahari terbenam lokal
                astrometric = observer_at.at(sunset_time_utc).observe(moon)
                astrometric_geocentric = observer_at_geocenter.at(sunset_time_utc).observe(moon)

                apparent = astrometric.apparent()
                # apparent_geocentric = astrometric_geocentric.apparent() # Tidak digunakan

                alt, az, _ = apparent.altaz()

                astrometric_sun = observer_at.at(sunset_time_utc).observe(sun)
                astrometric_sun_geocentric = observer_at_geocenter.at(sunset_time_utc).observe(sun)

                # elongasi geosentris
                elongation = astrometric_geocentric.separation_from(astrometric_sun_geocentric).degrees
                topo_moon = observer_at.at(sunset_time_utc).observe(moon).apparent()
                topo_sun = observer_at.at(sunset_time_utc).observe(sun).apparent()
                el_topo = topo_sun.separation_from(topo_moon).degrees

                # Simpan hasil
                altitudes[i, j] = alt.degrees
                elongations[i, j] = el_topo

    # Mengembalikan data yang dibutuhkan untuk plotting
    return {
        "latitudes": latitudes,
        "longitudes": longitudes,
        "altitudes": altitudes,
        "elongations": elongations,
        "jd_map": dates # Mengembalikan objek dates untuk digunakan di judul peta
    }

def plot_hilal_map(data_peta, bulan_hijriah_str, tahun_hijriah, data_type="altitude"):

    latitudes = data_peta["latitudes"]
    longitudes = data_peta["longitudes"]
    jd_map = data_peta["jd_map"]

    jd_caldat_result = caldat(float(jd_map.tt), 7).result

    # Tentukan data yang akan diplot, batas kontur, dan judul berdasarkan data_type
    if data_type.lower() == "altitude":
        plot_data = data_peta["altitudes"]
        bounds = [0, 3, 12]  # Nilai batas untuk ketinggian (sesuai kriteria)
        title_type = "KETINGGIAN HILAL (TOPOSENTRIK)"
    elif data_type.lower() == "elongation":
        plot_data = data_peta["elongations"]
        bounds = [0, 6.4, 15]  # Nilai batas untuk elongasi (sesuai kriteria NEO MABIMS)
        title_type = "ELONGASI (TOPOSENTRIK)"
    else:
        print(f"Error: data_type '{data_type}' tidak valid. Gunakan 'altitude' atau 'elongation'.")
        return # Keluar dari fungsi jika data_type tidak valid

    # --- Bagian Plotting ---

    fig = plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([94, 142, -12, 7], crs=ccrs.PlateCarree())

    # Add coastlines and borders with specified colors
    ax.coastlines('10m')

    ax.set_facecolor('white')
    ax.add_feature(cfeature.LAND, color='#BFD8A6')
    ax.add_feature(cfeature.OCEAN, color='white')
    ax.add_feature(cfeature.COASTLINE, edgecolor='#5A4D3D')
    ax.add_feature(cfeature.BORDERS, linestyle=':', edgecolor='#5A4D3D')
    ax.add_feature(cfeature.LAKES, color='#A2C4E6', alpha=0.6)
    ax.add_feature(cfeature.RIVERS, color='#76A5AF', alpha=0.7)

    # Tambahkan garis lintang dan bujur (gridlines) with specified styles
    gl = ax.gridlines(draw_labels=True, linestyle="--", linewidth=0.5, color="gray", alpha=0.7)
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 12, 'color': 'black'}
    gl.ylabel_style = {'size': 12, 'color': 'black'}

    # Contour plot using the selected data and bounds
    cmap = mcolors.ListedColormap(["red", "black"])
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    contour = ax.contour(longitudes, latitudes, plot_data, levels=7, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())

    # Add labels to the contour lines
    contour_labels = ax.clabel(contour, inline=True, fontsize=10, fmt='%.1f')

    # Add a white background to the labels for better visibility
    for label in contour_labels:
        label.set_backgroundcolor('white')


    ax.text(0.01, 0.03, "created by Institut Latihan Islam Malaysia (ILIM)", horizontalalignment='left', fontname="sans-serif",
                            verticalalignment='center', transform = ax.transAxes, fontsize=10, color='black')

    # Use dynamic title based on the data_type
    plt.title(f'PETA {title_type} WAKTU MATAHARI TERBENAM \n PENENTUAN AWAL BULAN {(hijriah().bulan_hijriah(bulan_hijriah_str)).upper()} {tahun_hijriah} H \n {jd_caldat_result.upper()} \n')
    return fig


    def __init__(self, bulan, tahun, lok, lat, lon, TZ='Asia/Kuala_Lumpur', TT=0, TH=0, kriteria='NEO MABIMS'):
        # Melakukan perhitungan awal bulan di dalam konstruktor
        self.awal_bulan_data = awalbulan(bulan, tahun, lok, lat, lon, TZ, TT, TH, kriteria).calculate_hilal()

        self.katsner = visibilitas_katsner(self.awal_bulan_data, lat, lon, TZ, TT)

    def tampilkan_data_hisab(self):
        """Menampilkan data hisab awal bulan."""
        printdata(self.awal_bulan_data)

    def tampilkan_data_katsner(self):
        return self.katsner

    def tampilkan_peta_tinggi_hilal(self):
        """Menampilkan peta tinggi hilal."""
        data_peta = datapeta(self.awal_bulan_data)
        plot_hilal_map(data_peta, bulan, tahun, data_type="altitude")

    def tampilkan_peta_elongasi_hilal(self):
        """Menampilkan peta elongasi hilal."""
        data_peta = datapeta(self.awal_bulan_data)
        plot_hilal_map(data_peta, bulan, tahun, data_type="elongation")

    def tampilkan_kurva_katsner(self, alat='TELESCOPE'):
        if self.katsner.empty:
            return
        else:
            plot_katsner_curve(self.katsner, self.awal_bulan_data, bulan, tahun, alat=alat)


st.markdown("<h1 style='text-align: center;'>HISAB AWAL BULAN HIJRIAH</h1>", unsafe_allow_html=True)

st.sidebar.header("Data Input")

# Input dari pengguna di sidebar
bulan_input = int(st.sidebar.number_input("Bulan Hijriah", min_value=1, max_value=12, value=6, step=1))
tahun_input = int(st.sidebar.number_input("Tahun Hijriah", min_value=1410, value=1447, step=1))
lokasi_input = st.sidebar.text_input("Nama Lokasi", "Kuala Lumpur")
lat_input = float(st.sidebar.number_input("Lintang (Decimal Degrees)", value=3.1388358,step=0.000001,format="%.6f"))
lon_input = float(st.sidebar.number_input("Bujur (Decimal Degrees)", value=101.522171,step=0.000001,format="%.6f"))
TZ_input = st.sidebar.selectbox("Timezone", ["Asia/Kuala_Lumpur"], index=0)
TT_input = int(st.sidebar.number_input("Ketinggian (m)", value=109,step=1))
TH_input = int(st.sidebar.number_input("Tambah Hari", value=0))
kriteria_input = st.sidebar.selectbox("Kriteria Kebolehnampakan", ["NEO MABIMS", "MABIMS LAMA"], index=0)
alat_katsner_input = st.sidebar.selectbox("Plot Model Kastner", ["NAKED EYE", "THEODOLIT", "TELESCOPE"], index=0)

awal_bulan_data = awalbulan(
    bulan_input, tahun_input, lokasi_input,
    lat_input, lon_input, TZ_input,
    TT_input, TH_input, kriteria_input
).calculate_hilal()

data = variabel(awal_bulan_data)

data_peta = datapeta(awal_bulan_data)
plot_altitude =  plot_hilal_map(data_peta, bulan_input, tahun_input, data_type="altitude")
plot_elongation =  plot_hilal_map(data_peta, bulan_input, tahun_input, data_type="elongation")
katsner = visibilitas_katsner(awal_bulan_data, lat_input, lon_input, TZ_input, TT_input)
if katsner.empty:
    pass
else:
    plot_katsner = plot_katsner_curve(katsner, awal_bulan_data,bulan_input, tahun_input, alat_katsner_input)
if st.sidebar.button("Hitung"):

    st.code(f"""
        --------------------------------------------------------------------------------
                                     DATA FALAK HILAL
                               AWAL BULAN {(hijriah().bulan_hijriah(bulan_input)).upper()} {tahun_input} H
                 JET PROPULSION LABORATORY (JPL) NASA, BY INSTITUT LATIHAN ISLAM MALAYSIA (ILIM)
        --------------------------------------------------------------------------------

        Pengiraan dilakukan untuk waktu matahari terbenam pada {data['sunset'].strftime('%H:%M:%S')}, tarikh {data['sunset'].day} {data['n_bln']} {data['sunset'].year} M
        Semua data dipaparkan dalam waktu tempatan pencerap
        Lokasi                   : {lokasi_input if lokasi_input else 'Tidak disebutkan'}
        Lintang                  : {konversi(lat_input, 'LINTANG').result}, Bujur: {konversi(lon_input, 'BUJUR').result}
        Elevasi                  : {TT_input:.2f} m
        Time Zone                : {TZ_input} {'+' if data['delta_time'] >= 0 else ''}{data['delta_time']}

        Visibilitas Oddeh
           Best Time             : {konversi(data['best_time'], 'JAM').result}
           q                     : {round(data['q'], 2)}
           Parameter             : {data['parameter']}

        Waktu Konjungsi          : {data['konjungsi'].day} {data['n1_bln']} {data['konjungsi'].year} M {data['konjungsi'].hour:02d}:{data['konjungsi'].minute:02d}:{data['konjungsi'].second:02d} LT
        Waktu Matahari Terbenam  : {data['sunset'].strftime('%H:%M:%S')}
        Waktu Bulan Terbenam     : {data['moonset'].strftime('%H:%M:%S')}
        Ketinggian Matahari      : {konversi(data['sun_alt']).result}
        Ketinggian Bulan         : {konversi(data['moon_alt']).result}
        Azimut Matahari          : {konversi(data['sun_az']).result}
        Azimut Bulan             : {konversi(data['moon_az']).result}
        Lebar Sabit Bulan        : {konversi(data['W']).result}
        Elongasi Bulan (Geo)     : {konversi(data['el_geo']).result}
        Elongasi Bulan (Topo)    : {konversi(data['el_topo']).result}
        Iluminasi Bulan          : {round(float(data['illumination']), 2)} %
        Umur Bulan               : {konversi(data['moonage'], 'JAM1').result}
        Lama Hilal di Ufuk       : {konversi(data['lag_time'], 'JAM1').result}
    """)

    st.markdown("<h1 style='text-align: center;'>PETA HILAL</h1>", unsafe_allow_html=True)

    st.pyplot(plot_altitude)
    st.pyplot(plot_elongation)

    if katsner.empty:
        pass 
    else:
        st.markdown("<h1 style='text-align: center;'>KELUK KEBOLEHNAMPAKAN MODEL KATSNER</h1>", unsafe_allow_html=True)
        st.pyplot(plot_katsner)



