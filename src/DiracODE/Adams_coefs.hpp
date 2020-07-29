#pragma once
#include <array>

//******************************************************************************
template <int N> struct AdamsCoefs {};
//------------------------------------------------------------------------------
template <> struct AdamsCoefs<8> {
  // Adams coefs:
  static constexpr int N = 8;
  static constexpr std::array<int, N> AMa = {
      -33953, 312874, -1291214, 3146338, -5033120, 5595358, -4604594, 4467094};
  static constexpr double AMd = 1.0 / 3628800;
  static constexpr int AMaa = 1070017;
  // Outward coefs:
  static constexpr int OIe[N][N] = {
      {-1338, 2940, -2940, 2450, -1470, 588, -140, 15},
      {-240, -798, 1680, -1050, 560, -210, 48, -5},
      {60, -420, -378, 1050, -420, 140, -30, 3},
      {-32, 168, -672, 0, 672, -168, 32, -3},
      {30, -140, 420, -1050, 378, 420, -60, 5},
      {-48, 210, -560, 1050, -1680, 798, 240, -15},
      {140, -588, 1470, -2450, 2940, -2940, 1338, 105},
      {-960, 3920, -9408, 14700, -15680, 11760, -6720, 2283}};
  static constexpr std::array<int, N> OIa = {-105, 15, -5, 3, -3, 5, -15, 105};
  static constexpr int OId = 840;
};
//------------------------------------------------------------------------------
template <> struct AdamsCoefs<7> {
  // Adams coefs:
  static constexpr int N = 7;
  static constexpr std::array<int, N> AMa = {1375,   -11351,  41499, -88547,
                                             123133, -121797, 139849};
  static constexpr double AMd = 1.0 / 120960;
  static constexpr int AMaa = 36799;
  // Outward coefs:
  static constexpr int OIe[N][N] = {
      {-609, 1260, -1050, 700, -315, 84, -10},
      {-140, -329, 700, -350, 140, -35, 4},
      {42, -252, -105, 420, -126, 28, -3},
      {-28, 126, -420, 105, 252, -42, 4},
      {35, -140, 350, -700, 329, 140, -10},
      {-84, 315, -700, 1050, -1260, 609, 60},
      {490, -1764, 3675, -4900, 4410, -2940, 1089}};
  static constexpr std::array<int, N> OIa = {-60, 10, -4, 3, -4, 10, -60};
  static constexpr int OId = 420;
};
//------------------------------------------------------------------------------
template <> struct AdamsCoefs<6> {
  // Adams coefs:
  static constexpr int N = 6;
  static constexpr std::array<int, N> AMa = {-863,  6312,   -20211,
                                             37504, -46461, 65112};
  static constexpr double AMd = 1.0 / 60480;
  static constexpr int AMaa = 19087;
  // Outward coefs:
  static constexpr int OIe[N][N] = {
      {-77, 150, -100, 50, -15, 2}, {-24, -35, 80, -30, 8, -1},
      {9, -45, 0, 45, -9, 1},       {-8, 30, -80, 35, 24, -2},
      {15, -50, 100, -150, 77, 10}, {-72, 225, -400, 450, -360, 147}};
  static constexpr std::array<int, N> OIa = {-10, 2, -1, 1, -2, 10};
  static constexpr int OId = 60;
};
//------------------------------------------------------------------------------
template <> struct AdamsCoefs<5> {
  // Adams coefs:
  static constexpr int N = 5;
  static constexpr std::array<int, N> AMa = {27, -173, 482, -798, 1427};
  static constexpr double AMd = 1.0 / 1440;
  static constexpr int AMaa = 475;
  // Outward coefs:
  static constexpr int OIe[N][N] = {{-65, 120, -60, 20, -3},
                                    {-30, -20, 60, -15, 2},
                                    {15, -60, 20, 30, -3},
                                    {-20, 60, -120, 65, 12},
                                    {75, -200, 300, -300, 137}};
  static constexpr std::array<int, N> OIa = {-12, 3, -2, 3, -12};
  static constexpr int OId = 60;
};
