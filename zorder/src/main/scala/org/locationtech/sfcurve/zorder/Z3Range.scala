/***********************************************************************
* Copyright (c) 2013-2015 Commonwealth Computer Research, Inc.
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Apache License, Version 2.0 which
* accompanies this distribution and is available at
* http://www.opensource.org/licenses/apache2.0.php.
*************************************************************************/
package org.locationtech.sfcurve.zorder

/**
 * Represents a cube in index space defined by min and max as two opposing points.
 * All operations refer to index space.
 */
case class Z3Range(min: Long, max: Long) extends ZRange {

  // contains in user space - each dimension is contained
  override def containsInUserSpace(bits: Long) = {
    val (x, y, z) = Z3(bits).decode
    x >= Z3(min).d0 && x <= Z3(max).d0 && y >= Z3(min).d1 && y <= Z3(max).d1 && z >= Z3(min).d2 && z <= Z3(max).d2
  }

  // overlap in user space - if any dimension overlaps
  override def overlapsInUserSpace(r: ZRange): Boolean =
    overlaps(Z3(min).d0, Z3(max).d0, Z3(r.min).d0, Z3(r.max).d0) &&
        overlaps(Z3(min).d1, Z3(max).d1, Z3(r.min).d1, Z3(r.max).d1) &&
        overlaps(Z3(min).d2, Z3(max).d2, Z3(r.min).d2, Z3(r.max).d2)

  private def overlaps(a1: Int, a2: Int, b1: Int, b2: Int) = math.max(a1, b1) <= math.min(a2, b2)
}
