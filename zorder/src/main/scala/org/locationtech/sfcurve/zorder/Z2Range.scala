/***********************************************************************
 * Copyright (c) 2015 Azavea.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Apache License, Version 2.0 which
 * accompanies this distribution and is available at
 * http://www.opensource.org/licenses/apache2.0.php.
 ***********************************************************************/

package org.locationtech.sfcurve.zorder

/**
 * Represents a rectangle in defined by min and max as two opposing points
 * 
 * @param min: lower-left point 
 * @param max: upper-right point
 */
case class Z2Range(min: Long, max: Long) extends ZRange {

  // contains in user space - each dimension is contained
  override def containsInUserSpace(bits: Long): Boolean = {
    val (x, y) = Z2(bits).decode
    x >= Z2(min).d0 && x <= Z2(max).d0 && y >= Z2(min).d1 && y <= Z2(max).d1
  }

  // overlap in user space - if any dimension overlaps
  override def overlapsInUserSpace(r: ZRange): Boolean =
    overlaps(Z2(min).d0, Z2(max).d0, Z2(r.min).d0, Z2(r.max).d0) &&
        overlaps(Z2(min).d1, Z2(max).d1, Z2(r.min).d1, Z2(r.max).d1)

  private def overlaps(a1: Int, a2: Int, b1: Int, b2: Int) = math.max(a1, b1) <= math.min(a2, b2)
}
