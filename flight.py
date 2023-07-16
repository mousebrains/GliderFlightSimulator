#! /usr/bin/env python3
#
# Fly a glider through Javier's ROMS velocity fields
#
# Since I'm only interested in surfacing positions and all the subsurface
# drifts are linear in time, step surfacing to surfacing
#
# July-2023, Pat Welch, pat@mousebrains.com

import xarray as xr
import pandas as pd
import numpy as np
from TPWUtils.GreatCircle import greatCircle as distance
import time
import sys

class NAV: # Fly a line of longitude
    def __init__(self, lon0:float, lat0:float, lon1:float, speed:float) -> None:
        self.__lon0 = lon0 # Initial longitude (deg)
        self.__lat0 = lat0 # Initial latitude (deg)
        self.speed = speed # Horizontal through water speed of glider (m/s)
        self.__xMax = float(distance(lon0, lat0, lon1, lat0)) # Meters from lon0 to lon1 along lat0
        self.__dLondDeg = self.__xMax / abs(lon1 - lon0) # m/deg of longitude
        lonMid = (lon0 + lon1) / 2
        self.__dLatdDeg = distance(lonMid, lat0-0.5, lonMid, lat0+0.5) # m/deg of latitude
        self.qWestward = lon1 < lon0
        if self.qWestward: self.__xMax = -self.__xMax

    def qDone(self, x, y) -> bool:
        return (self.qWestward and (x < self.__xMax)) or (not self.qWestward and (x > self.__xMax))

    def lon(self, x:float) -> float: return float(self.__lon0 + x / self.__dLondDeg)
    def lat(self, y:float) -> float: return float(self.__lat0 + y / self.__dLatdDeg)

    def uGlider0(self) -> float: return -self.speed if self.qWestward else self.speed
    def vGlider0(self) -> float: return 0

    def theta(self, x, y) -> float: return np.arctan2(y, self.__xMax - x)


class ConstantHeading(NAV):
    # Fly at a constant heading
    def __init__(self, lon0:float, lat0:float, lon1:float, speed:float, heading:float) -> None:
        NAV.__init__(self, lon0, lat0, lon1, speed)
        self.__velocity = (
                np.sin(np.deg2rad(heading)) * speed,
                np.cos(np.deg2rad(heading)) * speed,
                )

    def gliderVelocity(self, x:float, y:float, uDA:float, vDA:float, tTotal:float) -> tuple:
        return self.__velocity

    def uGlider0(self) -> float: return self.__velocity[0]
    def vGlider0(self) -> float: return self.__velocity[1]

class CurrentCorrection(NAV):
    # Current correction aiming for far end of line
    def __init__(self, lon0:float, lat0:float, lon1:float, speed:float) -> None:
        NAV.__init__(self, lon0, lat0, lon1, speed)
        self.__forward = (-speed if self.qWestward else speed, 0) # Full speed towards lon1

    def gliderVelocity(self, x:float, y:float, uDA:float, vDA:float, tTotal:float) -> tuple:
        theta = self.theta(x, y) # Angle in cartesian plan to (xMax,0)
        ctheta = np.cos(theta)
        stheta = np.sin(theta)
        uPrime = uDA * ctheta - vDA * stheta # Rotate Depth Averaged current along line to (xMax,0)
        vPrime = uDA * stheta + vDA * ctheta

        if abs(vPrime) > self.speed: return self.__forward # We can't do a ferry, to much current

        # We can ferry
        vGldPrime = -vPrime # Cancel out vPrime
        phi = np.arcsin(vGldPrime / self.speed) 
        uGldPrime = self.speed * np.cos(phi)
        uGld =  uGldPrime * ctheta + vGldPrime * stheta
        vGld = -uGldPrime * stheta + vGldPrime * ctheta
        return (uGld, vGld)

class FlyToLine(NAV):
    # Try and get back on the line "soonest"
    def __init__(self, lon0:float, lat0:float, lon1:float, speed:float) -> None:
        NAV.__init__(self, lon0, lat0, lon1, speed)
        self.__uGldMin = 0.1 * speed
        if self.qWestward: self.__uGldMin = -self.__uGldMin
        self.__vGldMax = np.sqrt(speed*speed - self.__uGldMin)

    def gliderVelocity(self, x:float, y:float, uDA:float, vDA:float, tTotal:float) -> tuple:
        vGld = -y / tTotal - vDA
        if abs(vGld) > self.speed:
            uGld = self.__uGldMin
            vGld = self.__vGldMax
            if y > 0: vGld = -vGld
            return (uGld, vGld)
        uGld = np.sqrt(self.speed*self.speed - vGld * vGld)
        if self.qWestward: uGld = -uGld
        return (uGld, vGld)

def flyLine(fn:str, stime:str, args:dict) -> pd.DataFrame:
    lonEast = args["lonEast"]
    lonWest = args["lonWest"]
    qWestward = args["qWestward"]
    dzdt = args["dzdt"]
    pitch = args["pitch"]
    aoa = args["aoa"]
    depthMax = args["depthMax"]
    nYos = args["nYos"]
    nIterations = args["iterationLimit"]
    dtSurfaceInt = args["dtSurface"]
    dtSurface = np.timedelta64(dtSurfaceInt, "s")

    spd = dzdt / np.sin(np.deg2rad(pitch + aoa)) # Glider's horizontal speed through water

    with xr.open_dataset(fn) as ds:
        stime = min(ds.t.data) if stime is None else np.datetime64(stime)

        if lonEast is None: lonEast = min(ds.lon.data)
        if lonWest is None: lonWest = max(ds.lon.data)

        depths = ds.depth.data
        timeToDepth = np.concatenate((depths, [depthMax])) / dzdt # dive time to depths
        dTimeToDepth = np.diff(timeToDepth) # dive time between depths

        iLonMin = min(ds.lon.data) # Data constraints for interpolation
        iLonMax = max(ds.lon.data)
        iTimeMin = min(ds.t.data)
        iTimeMax = max(ds.t.data)

        t = stime;
        lat = ds.lat.data[0] # The line is constant latitude
        lon0 = lonWest if qWestward else lonEast
        lon1 = lonEast if qWestward else lonWest

        if args["heading"] is not None: # Fly a constant heading
            nav = ConstantHeading(lon0, lat, lon1, spd, args["heading"])
        elif args["toLine"]: # Fly back to the line as soon as possible
            nav = FlyToLine(lon0, lat, lon1, spd)
        else: # Constant current is default
            nav = CurrentCorrection(lon0, lat, lon1, spd)

        lon = lon0 # Starting longitude

        # Starting point

        steps = [pd.DataFrame(dict(
            t=[pd.Timestamp(t)], days=0,
            lon=lon, lat=lat,
            x=0, y=0,
            uGld=nav.uGlider0(), vGld=nav.vGlider0(), # Initial glider velocities
            waterVx=0, waterVy=0,
            bottom=None, theta=0))]

        for iteration in range(nIterations): # Infinite loop prevention
            prev = steps[-1].iloc[-1]
            lon = nav.lon(prev.x)

            pt = dict(depth = depths,
                      t = min(max(t, iTimeMin), iTimeMax),
                      lon = min(max(lon, iLonMin), iLonMax))
            u = ds.u.interp(pt, assume_sorted=True).data # Water velocities at depths
            v = ds.v.interp(pt, assume_sorted=True).data

            u = np.concatenate([u, [u[-1]]]) # Repeat last value
            v = np.concatenate([v, [u[-1]]]) # Repeat last value

            uMid = (u[:-1] + u[1:]) / 2 # mean velocity in bin, assuming linear transit over bin
            vMid = (v[:-1] + v[1:]) / 2

            qMid = np.isnan(uMid) # below bottom, velocities are nan
            if np.any(qMid): # Bottom seen
                if np.all(qMid):
                    print("No valid depths found at", lon)
                    break
                uMid[qMid] = 0 # No current below the bottom
                vMid[qMid] = 0
                bottom = depths[np.logical_not(np.isnan(u[:-1]))][-1]
            else: # full water column depths
                bottom = depthMax

            tPerDive = bottom / dzdt # Time to dive
            tTotal = int(np.round(tPerDive * 2 * nYos)) # Total time between surfacing
            uDA = np.sum(dTimeToDepth * uMid) / tPerDive # Depth averaged current
            vDA = np.sum(dTimeToDepth * vMid) / tPerDive

            x = float(prev.x + (prev.uGld + uDA) * tTotal) # Where the glider surfaces
            y = float(prev.y + (prev.vGld + vDA) * tTotal)
            t += np.timedelta64(tTotal, "s") # Surfacing time

            (uGld, vGld) = nav.gliderVelocity(x, y, uDA, vDA, tTotal)

            dive = pd.Series(dict(t=pd.Timestamp(t),
                                  days=(t - stime).astype("timedelta64[s]").astype(float)/86400,
                                  lon=nav.lon(x), lat=nav.lat(y),
                                  x=x, y=y,
                                  uGld=uGld, vGld=vGld,
                                  waterVx=uDA, waterVy=vDA, bottom=bottom,
                                  theta=nav.theta(x, y)))
            # Surface drift
            t += dtSurface
            x += u[0] * dtSurfaceInt
            y += v[0] * dtSurfaceInt

            surf = pd.Series(dict(t=pd.Timestamp(t),
                                  days=(t - stime).astype("timedelta64[s]").astype(float)/86400,
                                  lon=nav.lon(x), lat=nav.lat(y),
                                  x=x, y=y,
                                  uGld=dive.uGld, vGld=dive.vGld,
                                  waterVx=u[0], waterVy=v[0], bottom=0,
                                  theta=nav.theta(x, y)))
            df = pd.DataFrame([dive, surf])
            steps.append(df)

            if nav.qDone(x, y): break # After surfacing

        # df = pd.concat(steps, ignore_index=True)
        # print(df)
        # sys.exit(1)
        return pd.concat(steps, ignore_index=True)

def flyDates(ifn:str, ofn:str, stime:str, etime:str, args:dict) -> pd.DataFrame:
    if stime is not None: stime = np.datetime64(stime)
    if etime is not None: etime = np.datetime64(etime)

    if stime is None or etime is None:
        with xr.open_dataset(ifn) as ds:
            if stime is None: stime = min(ds.t.data)
            if etime is None: etime = max(ds.t.data)

    items = []
    dt = None
    while stime <= etime:
        print(stime, etime)
        t = time.time()
        args["stime"] = stime
        df = flyLine(ifn, stime, args)
        df["stime"] = np.array([stime] * df.t.size)
        items.append(df)
        dt = time.time() - t
        print("Flew", stime, len(items), (etime - stime).astype("timedelta64[D]"), 
              "n", items[-1].shape[0], "dt", dt, "max y", max(abs(df.y)))
        stime += np.timedelta64(1, "D")
        print(stime, etime)

    df = pd.concat(items, ignore_index=True)
    print(df)
    ds = df.to_xarray()

    for key in ("heading", "toLine", "current", "qWestward", "lonEast", "lonWest", 
                "dzdt", "pitch", "aoa", "depthMax", "nYos", "dtSurface"):
        ds[key] = args[key]

    for key  in ds:
        ds[key].encoding.update(dict(zlib=True, complevel=5))

    print(ds)
    ds.to_netcdf(ofn)

if __name__ == "__main__":
    from argparse import ArgumentParser
    
    parser = ArgumentParser()
    parser.add_argument("input", type=str, help="Input NetCDF file")
    parser.add_argument("output", type=str, help="Output NetCDF file")
    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument("--heading", type=float, help="Fly glider at this heading (deg north)")
    grp.add_argument("--toLine", action="store_true", help="Fly glider back to the line")
    grp.add_argument("--current", action="store_true", help="Fly glider with current correction")
    parser.add_argument("--stime", type=str, default="2013-01-01",
                        help="Start flying lines on this date")
    parser.add_argument("--etime", type=str, default="2022-12-01",
                        help="Ending flying lines on this date")
    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument("--qWestward", action="store_true", help="Fly east to west")
    grp.add_argument("--qEastward", action="store_true", help="Fly west to east")
    parser.add_argument("--lonEast", type=float, default=109.43, help="Eastern longitude")
    parser.add_argument("--lonWest", type=float, default=111., help="Western longitude")
    parser.add_argument("--dzdt", type=float, default=0.1, help="Glider vertical fall rate (m/s)")
    parser.add_argument("--pitch", type=float, default=26., help="Glider's pitch (deg)")
    parser.add_argument("--aoa", type=float, default=2.5, help="Glider's angle-of-attack (deg)")
    parser.add_argument("--depthMax", type=float, default=350., help="Maximum dive depth (m)")
    parser.add_argument("--nYos", type=int, default=3, help="Number of yos between surfacing")
    parser.add_argument("--dtSurface", type=int, default=600, help="Seconds at the surface")
    parser.add_argument("--iterationLimit", type=int, default=1000,
                        help="Maximum number of surfacings")
    args = parser.parse_args()

    flyDates(args.input, args.output, args.stime, args.etime, vars(args))
