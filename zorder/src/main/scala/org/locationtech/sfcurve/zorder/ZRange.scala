/*
 * Copyright (c) 2013-2016 Commonwealth Computer Research, Inc.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Apache License, Version 2.0 which
 * accompanies this distribution and is available at
 * http://www.opensource.org/licenses/apache2.0.php.
 */

package org.locationtech.sfcurve.zorder

trait ZRange[ZRangeT <: ZRange[ZRangeT]] {

  require(min <= max, s"Not: $min < $max")

  private [zorder] def zObject: ZObject[ZRangeT]

  def min: Long
  def max: Long

  def mid: Long = (min + max) >>> 1 // overflow safe mean

  def length: Long = max - min + 1

  // contains in index space (e.g. the long value)
  private [sfcurve] def untypedContains(z: Long): Boolean = z >= min && z <= max

  // contains in index space (e.g. the long value)
  def contains(r: ZRangeT): Boolean = untypedContains(r.min) && untypedContains(r.max)

  // overlap in index space (e.g. the long value)
  def overlaps(r: ZRangeT): Boolean = untypedContains(r.min) || untypedContains(r.max)

  // contains in user space - each dimension is contained
  private [sfcurve] def untypedContainsInUserSpace(z: Long): Boolean = {
    var dim = 0
    while (dim < zObject.Dimensions) {
      val zI = zObject.decode(z, dim)
      val minI = zObject.decode(min, dim)
      if (zI < minI) {
        return false
      }
      val maxI = zObject.decode(max, dim)
      if (zI > maxI) {
        return false
      }
      dim += 1
    }
    true
  }

  // contains in user space - each dimension is contained
  def containsInUserSpace(r: ZRangeT): Boolean =
    untypedContainsInUserSpace(r.min) && untypedContainsInUserSpace(r.max)

  // overlap in user space - if any dimension overlaps
  def overlapsInUserSpace(r: ZRangeT): Boolean = {
    var dim = 0
    while (dim < zObject.Dimensions) {
      val minI = zObject.decode(min, dim)
      val maxI = zObject.decode(max, dim)
      val rMin = zObject.decode(r.min, dim)
      val rMax = zObject.decode(r.max, dim)
      if (math.max(minI, rMin) > math.min(maxI, rMax)) {
        return false
      }
      dim += 1
    }
    true
  }
}

object ZRange {

  def apply(min: Z2, max: Z2) = Z2Range(min.z, max.z)
  def apply(min: Z3, max: Z3) = Z3Range(min.z, max.z)

  private [sfcurve] def apply(min: Long, max: Long, dimensions: Int) = {
    if (dimensions == 2) {
      Z2Range(min, max)
    } else if (dimensions == 3) {
      Z3Range(min, max)
    } else {
      throw new NotImplementedError(s"No ZRange defined for $dimensions dimensions")
    }
  }
}

case class Z2Range(min: Long, max: Long) extends ZRange[Z2Range] {

  override private [zorder] val zObject = Z2

  // Z2 specific methods - we can't generalize these without incurring allocation costs
  def zdivide(z: Z2): (Z2, Z2) = Z2.zdivide(z, Z2(min), Z2(max))
  def contains(z: Z2): Boolean = untypedContains(z.z)
  def containsInUserSpace(z: Z2): Boolean = untypedContainsInUserSpace(z.z)
}

case class Z3Range(min: Long, max: Long) extends ZRange[Z3Range] {

  override private [zorder] val zObject = Z3

  // Z3 specific methods - we can't generalize these without incurring allocation costs
  def zdivide(z: Z3): (Z3, Z3) = Z3.zdivide(z, Z3(min), Z3(max))
  def contains(z: Z3): Boolean = untypedContains(z.z)
  def containsInUserSpace(z: Z3): Boolean = untypedContainsInUserSpace(z.z)
}

